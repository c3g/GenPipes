#!/usr/bin/env python

import argparse
import csv
import http.client
import os
import sys


def main():
    parser = argparse.ArgumentParser(description='Parse WGS project')
    parser.add_argument('-d', '--dir', help='Directory to extract from', required=True)
    parser.add_argument('-s', '--sample', help='Sample name', required=True)
    parser.add_argument('-a', '--arrayReport', help='Array Sample report', required=True)
    parser.add_argument('-v', '--vertical', help='Vertival Output', required=False, type=bool)
    parser.add_argument('-w', '--host', help='Nanuq host to send data', required=False)
    parser.add_argument('-x', '--authFile', help='Authentification file for connecting to host', required=False)
    parser.add_argument('-m', '--manifest', help='Version of the Array Manifest (ex. 1.2)', required=False, default="NA" )

    args = parser.parse_args()

    if not os.path.isdir(args.dir):
      print(args.dir + " is not a directory")
      sys.exit(1)

    print("Reading " + args.dir)

    sampleStats = {'dupPct':'','chimeras':'','meanISize':'','iWidth99':'','interQ75_Q50Depth':'','depth':'','identity':"",'callRate':''}
    getDups(args.sample, args.dir, sampleStats)
    getBarcode(args.sample, args.dir, sampleStats)
    getChimeras(args.sample, args.dir, sampleStats)
    getInsertSize(args.sample, args.dir, sampleStats)
    getDepth(args.sample, args.dir, sampleStats)
    getArrayIdentity(args.sample, args.dir, args.manifest, sampleStats)
    getArrayCallRate(args.sample, args.arrayReport, sampleStats)

    if args.vertical:
      print("MeanInsertSize     : " + sampleStats['meanISize'])
      print("99% InsertSize     : " + sampleStats['iWidth99'])
      print("Chimeras           : " + sampleStats['chimeras'])
      print("Dups               : " + sampleStats['dupPct'])
      print("Depth              : " + sampleStats['depth'])
      print("InterQuantile Depth: " + sampleStats['interQ75_Q50Depth'])
      print("Identity           : " + sampleStats['identity'])
      print("Call Rate          : " + sampleStats['callRate'])
      print("Sample             : " + passFail(sampleStats))
    else:
      print("Chimera rate,%Dups,InsertSize,Bin Width that contains 99% of Fragments Insert Size,InterQuantile Coverage Depth,Array call rate,% identity,coverage,PASS FAIL")
      outLine = sampleStats['chimeras'] + "," + sampleStats['dupPct'] + "," + sampleStats['meanISize'] + "," + sampleStats['iWidth99'] + "," + sampleStats['interQ75_Q50Depth'] + "," + sampleStats['callRate'] + "," + sampleStats['identity'] + "," + sampleStats['depth'] + "," + passFail(sampleStats)
      print(outLine)

    if args.host:
        url = '/nanuqMPS/changeProcessingStateWS?barcode={barcode}&meanISize={meanISize}&iWidth99={iWidth99}' \
              '&chimeras={chimeras}&dupPct={dupPct}&depth={depth}&interQ75_Q50Depth={interQ75_Q50Depth}' \
              '&identity={identity}&callRate={callRate}&passFail={passFail}' \
            .format(
                barcode=sampleStats['barcode'],
                meanISize=sampleStats['meanISize'],
                iWidth99=sampleStats['iWidth99'],
                chimeras=sampleStats['chimeras'],
                dupPct=sampleStats['dupPct'],
                depth=sampleStats['depth'],
                interQ75_Q50Depth=sampleStats['interQ75_Q50Depth'],
                identity=sampleStats['identity'],
                callRate=sampleStats['callRate'],
                passFail=passFail(sampleStats)
        )
        sys.exit(contact_server(args.host, url, auth_file=args.authFile) != 200)

def passFail(sampleStats):
    retVal=""

    if float(sampleStats['depth']) < 30:
        retVal += "LowDepth "
    if float(sampleStats['identity']) < 0.90:
        retVal += "IdentityLOW "
    if float(sampleStats['meanISize']) < 50:
        retVal += "InsertSizeLOW "
    if float(sampleStats['meanISize']) > 2000:
        retVal += "InsertSizeHIGH "
    if float(sampleStats['chimeras']) > 0.1:
        retVal += "chimericHIGH "
    if float(sampleStats['callRate']) < 0.90 and float(sampleStats['callRate']) > 0.0 :
        retVal += "callRateLOW "
    elif float(sampleStats['callRate']) == 0.0 :
        retVal += "callRateMISSING "

    if len(retVal) == 0:
        return "OK"
    return retVal

def getArrayCallRate(sample, report, sampleStats):
    arrayMetrics = csv.DictReader(open(report, 'rb'), delimiter=',')
    sample_not_found = True
    for line in arrayMetrics:
        if line['Name'] == sample:
            sample_not_found = False
            sampleStats['callRate'] = str(float(line['Call Rate'].translate(None, "% "))/100)
            break
    if sample_not_found:
        print("---------\nWARNING sample ("+ sample + ") is not found in the array sample report\n------\n")
        sampleStats['callRate'] = str(0.0)


def getArrayIdentity(sample, dir, ver, sampleStats):
    version= "_v" + ver if ver != "NA" else "" 
    snp_array_file='.snpArrayCmp'+ version +'.txt'
    MATCH_STR = '% Match : '
    with open(dir + '/' + sample + snp_array_file) as arrayMetrics:
        for line in arrayMetrics:
            if line.startswith(MATCH_STR):
                sampleStats['identity'] = str(float(line[len(MATCH_STR):].rstrip('\n'))/100)
                break

def getDepth(sample, dir, sampleStats):
    with open(dir + '/' + sample + '.sorted.dup.recal.coverage.tsv') as depthMetrics:
        depthMetrics.readline()
        line = depthMetrics.readline()
        values = line.split('\t')
        sampleStats['depth'] = values[8]
        sampleStats['interQ75_Q50Depth'] = str(int(values[9]) - int(values[11]))

def getInsertSize(sample, dir, sampleStats):
    with open(dir + '/' + sample + '.sorted.dup.recal.all.metrics.insert_size_metrics') as iSizeMetrics:
        for line in iSizeMetrics:
            if line.startswith('MEDIAN_INSERT_SIZE'):
                break
        for line in iSizeMetrics:
            if line == '\n':
                break
            values = line.split('\t')
            if not sampleStats['iWidth99']:
                sampleStats['iWidth99'] = values[17]
            else:
                sampleStats['iWidth99'] += ',' + values[17]

            if not sampleStats['meanISize']:
                sampleStats['meanISize'] = values[4]
            else:
                sampleStats['meanISize'] += ',' + values[4]

def getChimeras(sample, dir, sampleStats):
    with open(dir + '/' + sample + '.sorted.dup.recal.all.metrics.alignment_summary_metrics') as alnMetrics:
        for line in alnMetrics:
            if line.startswith('PAIR'):
                values = line.split('\t')
                if not sampleStats['chimeras']:
                    sampleStats['chimeras'] = values[20]
                else:
                    sampleStats['chimeras'] += ',' + values[20]

def getDups(sample, dir, sampleStats):
    with open(dir + '/' + sample + '.sorted.dup.metrics') as dupMetric:
        for line in dupMetric:
            if line.startswith('LIBRARY'):
                break
        for line in dupMetric:
            if line == '\n':
                break
            values = line.split('\t')
            if not values[7] == '?':
                if not sampleStats['dupPct']:
                    sampleStats['dupPct'] = values[7]
                else:
                    sampleStats['dupPct'] += ',' + values[7]

def getBarcode(sample, dir, sampleStats):
    with open(dir + '/' + sample + '.sorted.dup.metrics') as dupMetric:
        for line in dupMetric:
            if line.startswith('LIBRARY'):
                break
        for line in dupMetric:
            if line == '\n':
                break
            values = line.split('\t')
            sampleStats['barcode'] = values[0]

def contact_server(host, url, auth_file=None):
    https_connection = http.client.HTTPSConnection(host)
    https_connection.set_debuglevel(0)
    headers = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain"}

    if auth_file:
        with open(auth_file) as auth_file_handle:
            https_connection.request("POST", url, auth_file_handle, headers)
    else:
        https_connection.request("POST", url, None, headers)
    http_response = https_connection.getresponse()

    https_connection.close()
    return http_response.status

if __name__ == "__main__":
    main()
