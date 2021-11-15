#!/usr/bin/python
import time
import collections
import json
import re
import cgi
import datetime

#Put the path of the source mugqic_pipelines.log file
PIPELINE_DATA = "./mugqic_pipelines_statistics.log"

Interval = collections.namedtuple('Interval', ["date_from", "date_to"])
UsageRecord = collections.namedtuple('UsageRecord', ["date", "agent", "hostname", "pipeline", "steps", "nb_samples"])


def main():
    arguments = cgi.FieldStorage()

    #Merge all version of the same pipeline if requested
    if 'mergePipelineVersions' in arguments and arguments['mergePipelineVersions'].value == '1':
        mergePipelineVersions = True
    else:
        mergePipelineVersions = False

    #Only keep records starting at this date
    if 'date_from' in arguments:
        dateFrom = arguments['date_from'].value
    else:
        dateFrom = None

    #Only keep records up to this date
    if 'date_to' in arguments:
        dateTo = arguments['date_to'].value
    else:
        dateTo = None

    print('Content-Type: application/json')
    print('')
    print(getStatisticsJson(dateFrom, dateTo, mergePipelineVersions))


def getStatisticsJson(dateFrom, dateTo, mergePipelineVersions):
    rndDigits = 1

    #Categorize hostname to an HPC
    hostnameToHpcRegexp = {
        "^abacus*": "Abacus",
        "^lg-*": "Guillimin",
        "^lm-*": "Guillimin",
        "^ip\\d\\d": "Mammouth",
        "^briaree*": "Briaree"
    }
    interval = setInterval(dateFrom, dateTo)

    totalNbSamples = 0
    totalJobSubmissions = 0
    hostDict = {}
    pipelineDict = {}
    hostPipelineDict = {}
    earliestParsedDate = interval.date_to  # Store the date of the earliest record in the log file
    latestParsedDate = interval.date_from  # Store the date of the most recent record in the log file

    logfile = open(PIPELINE_DATA, "r")

    # Calculate global statistics, and split records per hostname, pipeline, or both at once
    for line in logfile:
        tokens = line.split("\t")
        timeStr = tokens[0]

        try:
            parsedDate = time.strptime(timeStr[:10], "%Y-%m-%d")

            if interval.date_from <= parsedDate <= interval.date_to:
                data_dict = {"date": parsedDate}
                for token in tokens[1:]:
                    keyVal = token.split("=")
                    data_dict[keyVal[0]] = keyVal[1].rstrip('\n')
                if data_dict["nb_samples"] == "":
                    continue

                # Some pipeline versions use different cases for the same pipeline.
                # dict["pipeline"] = dict["pipeline"].lower()

                # Merge all versions of a pipeline if requested, by stripping everything following an "-"
                if mergePipelineVersions:
                    data_dict["pipeline"] = data_dict["pipeline"].split('-', 1)[0]

                totalNbSamples += int(data_dict["nb_samples"])
                totalJobSubmissions += 1

                hostDict.setdefault(data_dict["hostname"], []).append(data_dict)
                pipelineDict.setdefault(data_dict["pipeline"], []).append(data_dict)
                hostPipelineDict.setdefault(data_dict["hostname"] + "|" + data_dict["pipeline"], []).append(data_dict)

                if parsedDate < earliestParsedDate:
                    earliestParsedDate = parsedDate
                if parsedDate > latestParsedDate:
                    latestParsedDate = parsedDate
        except:
            print("Warning: Pipelines log file has lines with invalid format:", line)
            continue

    # Prepare a hostname to HPC dictionary
    hostnameToHpc = {}
    for key in hostDict.keys():
        hostnameToHpc[key] = ""

        for testedRe in hostnameToHpcRegexp.keys():
            isMatch = re.match(testedRe, key)
            if isMatch:
                hpcName = hostnameToHpcRegexp[testedRe]
                hostnameToHpc[key] = hpcName
                break

    # Process records per hostname
    countPerHostDict = []
    for key in sorted(hostDict.keys(), key=lambda s: s.lower()):
        nbSamples = 0
        jobSubmissions = 0

        for record in hostDict[key]:
            nbSamples += int(record["nb_samples"])
            jobSubmissions += 1

        resultDict = {
            "hostname": key,
            "hpc": hostnameToHpc[key],
            "nb_submissions": jobSubmissions,
            "nb_samples": nbSamples,
            "average_samples_per_job": round(nbSamples / jobSubmissions, rndDigits)
        }

        countPerHostDict.append(resultDict)

    # Process records per pipeline
    countPerPipelineDict = []
    for key in sorted(pipelineDict.keys(), key=lambda s: s.lower()):
        nbSamples = 0
        jobSubmissions = 0

        for record in pipelineDict[key]:
            nbSamples += int(record["nb_samples"])
            jobSubmissions += 1

        resultDict = {
            "pipeline_name": key,
            "nb_submissions": jobSubmissions,
            "nb_samples": nbSamples,
            "average_samples_per_job": round(nbSamples / jobSubmissions, rndDigits)
        }

        countPerPipelineDict.append(resultDict)

    # Process records per hostname/pipeline
    countPerHostPipelineDict = []
    for key in sorted(hostPipelineDict.keys(), key=lambda s: s.lower()):
        nbSamples = 0
        jobSubmissions = 0

        hostname = hostPipelineDict[key][0]["hostname"]
        pipeline = hostPipelineDict[key][0]["pipeline"]

        for record in hostPipelineDict[key]:
            nbSamples += int(record["nb_samples"])
            jobSubmissions += 1

        resultDict = {
            "hostname": hostname,
            "hpc": hostnameToHpc[hostname],
            "pipeline_name": pipeline,
            "nb_submissions": jobSubmissions,
            "nb_samples": nbSamples,
            "average_samples_per_job": round(nbSamples / jobSubmissions, rndDigits)
        }

        countPerHostPipelineDict.append(resultDict)

    if totalJobSubmissions > 0:
        avgSamplesPerJob = round(totalNbSamples / totalJobSubmissions, rndDigits)
        startString = time.strftime("%Y-%m-%d", earliestParsedDate)
        endString = time.strftime("%Y-%m-%d", latestParsedDate)
    else:
        avgSamplesPerJob = 0
        startString = time.strftime("%Y-%m-%d", interval.date_from)
        endString = time.strftime("%Y-%m-%d", interval.date_to)

    # Build JSON hash
    jsonDict = {}
    jsonDict["settings"] = {
        "nb_submissions": totalJobSubmissions,
        "nb_samples": totalNbSamples,
        "average_samples_per_job": avgSamplesPerJob,
        "start": startString,
        "end": endString
    }
    jsonDict["per_host"] = countPerHostDict
    jsonDict["per_pipeline"] = countPerPipelineDict
    jsonDict["per_host_and_pipeline"] = countPerHostPipelineDict

    return json.dumps(jsonDict)


def setInterval(dateFromStr, dateToStr):
    fromDate = time.strptime("2014-01-01", "%Y-%m-%d")
    toDate = time.localtime()

    if dateFromStr is not None:
        try:
            fromDate = datetime.datetime.strptime(dateFromStr, "%Y-%m-%d").timetuple()
        except Exception as e:
            print e
            print("[Pipelines Usage Log] Invalid fromDate:", dateFromStr)
            pass

    if dateToStr is not None:
        try:
            toDate = datetime.datetime.strptime(dateToStr, "%Y-%m-%d").timetuple()
        except:
            print("[Pipelines Usage Log] Invalid toDate:", dateToStr)
            pass

    return Interval(date_from=fromDate, date_to=toDate)


if __name__ == "__main__":
    main()
