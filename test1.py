"""Usage: test.py (-h --help)  ...
       test.py -t <testType> ...
        
testType:
  1 variantCallingCore test
  2 variantCalling test
  3 addVarinatToIntVcArray test
Options:
  -h --help      | help
  -t testType    | execute specified test
"""

import variant_calling
import msvcrt
from docopt import docopt

def testVariantCallingCoreFunction():
    print("Executing variantCallingCore function test...")
    variant_calling.VCFfile = open(variant_calling.dir_path + "\output\VCF_test.vcf", "a+")
    variant_calling.VCFfile.seek(0)
    variant_calling.VCFfile.truncate()
    variant_calling.VCFfile.seek(0)
    variant_calling.content = [["1", "336", "C", "34", "A$A$A$A$AAAA,.,AA,,..AAAAaaaa...CcTtGg", "test"], ["1", "400", "A", "12", "...,,+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC+2AC,.,...,", "test"], ["1", "500", "A", "12", "...,,-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC-2AC,.,...,", "test"], ["1", "501", "A", "1", ".", "test"], ["1", "502", "C", "1", ".", "test"]]
    variant_calling.variantCalling()
    variant_calling.VCFfile.seek(0)
    #print(VCFfile.read())
    if variant_calling.VCFfile.read() == "1	336	C	A	GT	0/1\n1	400	A	AAC	GT	1/1\n1	500	AAC	A	GT	1/1\n":
        print("Test run succesufully...")
    else:
        print("Test fail...")
        
def testVariantCallingFunction():
    print("Executing variantCalling function test...")
    incVcArray = variant_calling.variantCallingCore(["1", "336", "C", "34", "A$A$A$A$AAAA,.,.,,,..AAAAaaaa...CcTtGg+2AC", "test"])
    sum = 0
    print(incVcArray)
    for x in incVcArray:
        sum += x[0]
    print(sum)
    if sum != 35:
        print("Test fail...")
    if incVcArray[0][0] == 12 and incVcArray[1][0] == 2 and incVcArray[2][0] == 2 and incVcArray[3][0] == 16 and incVcArray[4][0] == 2 and incVcArray[5][0] == 1:
        print("Test run succesufully...")
    else:
        print("Test fail...")
        
def testAddVarinatToIntVcArray():
    print("Executing addVarinatToIntVcArray function test...")
    p = ['.', 'T', 'G', 'A', 'C']
    incVcArray = [[0, p[i]] for i in range(len(p))]
    variant_calling.addVarinatToIntVcArray(incVcArray, "+2AC")
    if incVcArray[5][1] == "+2AC":
        print("Insertion added succesufully...")
        variant_calling.addVarinatToIntVcArray(incVcArray, "-2AT")
        if incVcArray[6][1] == "-2AT":
            print("Deletion added succesufully...")
            print("Function addVarinatToIntVcArray executed succesufully...")
        else:
            print("Test failed due to deletion adding!")
    else:
        print("Test failed due to insertion adding!")
            
def testFunction(m):

    global VCFfile  
            
    if str(m) == '1':
        testVariantCallingCoreFunction()

    elif str(m) == '2':
        testVariantCallingFunction()
        
    elif str(m) == '3':
        testAddVarinatToIntVcArray()

    else:
        print("Entered value is not allowed...")
        
# Run the program
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    testFunction(arguments['-t'][0])
