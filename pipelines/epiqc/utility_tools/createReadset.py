import os

"""
Searces for a directory named bigwig containing bigwig files in the working directory and creates
an readset from those files.
"""

def make_readset():
    file = open("readset.tsv", "w+")
    file.write("Sample\tReadset\tLibrary\tRunType\tRun\tLane\tAdapter1\tAdapter2\tQualityOffset\tBED\tFASTQ1\tFASTQ2\tBAM\tUMI\tBIGWIG\n")
    cpt = 1
    for line in os.listdir("bigwig"):
        split = line.split(".")
        #filename = line[5:]
        sample = split[2]
        file.write(sample+"\treadset"+str(cpt)+"\t\tSINGLE_END\t1\t1\t\t\t\t\t\t\t\t\tbigwig/"+line+"\n")
        cpt += 1
    
    file.close()


if __name__ == "__main__":
    make_readset()
