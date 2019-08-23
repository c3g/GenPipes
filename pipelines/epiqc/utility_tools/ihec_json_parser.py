#!/usr/bin/env python

import json
import os
import shutil
import ssl
import urllib2
import threading
import time
import argparse
import logging

def parseJson(file):
    with open(file) as json_file:
        data = json.load(json_file)

        marks = []
        bigwig_files = []
        sample_names = []

        for sample in data['datasets']:
            marks.append(data['datasets'][sample]['ihec_data_portal']['assay'])
            bigwig_files.append(data['datasets'][sample]['browser']['signal_unstranded'][0]['big_data_url'])
            sample_names.append(data['datasets'][sample]['sample_id'])

        assert len(sample_names) == len(marks) == len(bigwig_files)

        return [marks, bigwig_files, sample_names]


def downloadFiles(url, output_path):
    """
        donwload_path is the folder where files will be put
    """
    try:
        context = ssl._create_unverified_context()
        response = urllib2.urlopen(url, context=context)
    except urllib2.HTTPError as e:
        error_file = open("error_files.txt", "a")
        error_file.write("Couldn't download file for url : " + url + "\n")
        print("The server couldn't fulfill the request.")
        raise Exception("Error code: " + str(e.code) +" "+url)
    except urllib2.URLError as e:
        error_file = open("error_files.txt", "a")
        error_file.write("Couldn't download file for url : " + url + "\n")
        print("We failed to reach a server !")
        raise Exception("reason : " + str(e.reason) + url)

    filename = url.split("/")[-1]
    complete_file_name = os.path.join(output_path, filename)
    
    writeFile = open(complete_file_name, "w+")
    # log.info("READING " + filename)
    file = response.read()
    log.info("WRITING " + filename)
    writeFile.write(file)
    writeFile.close()
    log.info("DONE :" + filename)

def createReadsetMarkList(marks, bigwig_files, sample_names, dataset_name):
    bigwig_filenames = []
    for bigwig in bigwig_files:
        bigwig = bigwig.split("/")[-1]
        bigwig_filenames.append(os.path.join(dataset_name, bigwig))

    readset_file = open("readset.tsv", "w+")
    readset_file.write("Sample\tReadset\tLibrary\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM\tUMI\tBIGWIG\n")

    mark_list = open("mark_list.txt", "w+")
    mark_list.write("Copy paste the mark list below into the marks field of the chromimpute section in epiqc.base.ini :\n\n")
    mark_list_line = ""

    for i in range(len(marks)):
        mark_list_line += marks[i] + ","
        readset_file.write(sample_names[i]+"\treadset"+str(i+1)+"\t\tSINGLE_END\t1\t1\t\t\t\t\t\t\t\t\t"+bigwig_filenames[i]+"\n")

    mark_list.write(mark_list_line[:-1])
    mark_list.close()
    readset_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Downloads bigwig files from a json file then creates a readset and a string of marks to be used for epiqc")
    parser.add_argument("-j", "--jsonfile", help="File containing json information", type=str, required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to be used", type=int, required=False, default=4)

    args = parser.parse_args()

    log = logging.getLogger("Download data")
    log.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(format)
    log.addHandler(ch)

    json_output = parseJson(args.jsonfile)

    dataset_name = "bigwig"
    if os.path.exists(dataset_name):
        shutil.rmtree(dataset_name)
    os.mkdir(dataset_name)

    print("Starting bigwig download")
    start = time.time()
    threads = [threading.Thread(target=downloadFiles, args=(url.encode("utf-8"),dataset_name,)) for url in json_output[1]]

    NUMBER_OF_THREADS = args.threads #Downloads the files in batches threads                                                
    for i in range(0, len(threads), NUMBER_OF_THREADS):
            subthreads = threads[i:i+NUMBER_OF_THREADS]
            for thread in subthreads:
                    thread.start()
            for thread in subthreads:
                    thread.join()

    print("Elapsed time : " + str(time.time() - start))

    createReadsetMarkList(json_output[0], json_output[1], json_output[2], dataset_name)


