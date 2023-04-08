#!/usr/bin/env python3
import multiprocessing, time
from multiprocessing import Semaphore, Queue, Process, Array
from streaming.dump_reader import MultipleRecordReader


def startProcesses(processes):
    for process in processes:
        process.start()


def recordWorker(processId, recordkeeper, argsForEachRecord):
    while 1:
        # Wait for data to be available to be consumed
        recordkeeper.consumer_sem.acquire()

        # Get record from input queue
        record = recordkeeper.inputQueue.get()

        # Signal input reader that a spot in the input queue is free to be filled
        recordkeeper.producer_sem.release()

        # Check for poison pill
        if record is None:
            recordkeeper.worker_sem.release()
            break

        # Process records
        task_result = recordkeeper.record_handler(record, argsForEachRecord)

        # Construct result pair
        result_data = (record.get_timestep(), task_result)

        # Insert into output queue
        recordkeeper.outputQueue.put(result_data)


class MultipleRecordProcessor:
    # a record contains all of the info
    # relative to a single timestep; it is made of multiple entries

    def __init__(self, record_reader, data):
        # The object responsible for reading records from the input file
        # MultipleRecordReader(open(self.dumpFilePath, "r"), data)
        self.dumps_reader = record_reader

        # Semaphores for input queue synchronization
        self.consumer_sem = None
        self.producer_sem = None

        # Semaphore for signaling termination of worker threads
        self.worker_sem = None

        # Worker routine for each record
        self.record_handler = None

    def processEntriesInParallel(
            self, record_handler, task_additional_args, num_processes=0
    ):
        self.initialize(record_handler, task_additional_args, num_processes)
        processes = self.createProcesses(recordWorker)
        startProcesses(processes)
        print("Started processes.")
        print("Reading input into queue...")
        self.populateInputQueues()
        print("Input fully read.")
        print("Waiting for processing to finish")
        self.addPoisonPillForEachProcess()
        self.waitForProcessesToComplete()
        print("All processes have completed")

    def initialize(self, record_handler, task_additional_args, num_processes=0):
        self.initAttributes(record_handler, task_additional_args)
        self.initQueues()
        self.initNumProcesses(num_processes)
        self.initSynchronization()

    def initAttributes(self, record_handler, task_additional_args):
        self.record_handler = record_handler
        self.record_handler_args = [self, task_additional_args]

    def initQueues(self):
        self.inputQueue = Queue()
        self.outputQueue = Queue()

    def initNumProcesses(self, num_processes=0):
        if num_processes != 0:
            self.numProcesses = num_processes
        else:
            self.numProcesses = multiprocessing.cpu_count()

    def initSynchronization(self):
        # Reserve limited amount of queue spots for input data
        self.producer_sem = Semaphore(2 * self.numProcesses)

        # Signal for consumers to be notified of input availability
        self.consumer_sem = Semaphore(0)

        # Semaphore for signaling that a worker has terminated
        self.worker_sem = Semaphore(0)

    def initEndingArray(self):
        self.endingArray = Array("i", self.numProcesses)
        for i in range(self.numProcesses):
            self.endingArray[i] = 0

    def createProcesses(self, consumeRecord):
        processes = []
        for i in range(self.numProcesses):
            processes.append(
                Process(target=consumeRecord, args=(i, self, self.record_handler_args))
            )
        return processes

    # Successively read the input records to be processed
    # Will wait if queue is too full to prevent memory overflow
    def populateInputQueues(self):
        while 1:
            # Get record from record reader
            record = self.dumps_reader.readRecord()

            # Check if end of record keeping file is reached
            if record is None:
                break
            else:
                # Append record to queue
                self.insertRecordIntoQueue(record)

    # Method to insert records into worker queue if queue is not yet overflowing
    # Will wait if queue is full until a worker has processed at least one entry
    def insertRecordIntoQueue(self, queue_entry):
        # Ensure that there is enough space in the queue to insert an entry by reserving one spot
        self.producer_sem.acquire()

        # Append the data
        self.inputQueue.put(queue_entry)

        # Signal the consumers that an entry is available
        self.consumer_sem.release()

    def addPoisonPillForEachProcess(self):
        for i in range(self.numProcesses):
            print("Poison pill #{0} added".format(i))
            self.addPoisonPill()

    def addPoisonPill(self):
        self.insertRecordIntoQueue(None)

    def waitForProcessesToComplete(self):
        for i in range(self.numProcesses):
            self.worker_sem.acquire()
