import os

def main ():

    for errFile in os.popen("find batch_trees/condor_logs/*ttV_Jan23*.stderr -size +0k -print"):
        sampleLines = open(errFile.split('.',1)[0] + '.stdout').read().splitlines()
        if len(sampleLines) < 3:
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            print "echo 'No sample name in %s'" % (errFile.split('.',1)[0] + '.stdout')
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            continue
        if 'Passed cuts' in sampleLines[len(sampleLines)-2]:
            #print '%s exists, but .stdout file has "Passed cuts"' % errFile 
            continue

        index = 0
        if '.crc.nd.edu' in sampleLines[0]:
            index = 1
        sampleNameLine = sampleLines[2+index]
        if 'sample name =' in sampleNameLine:
            sampleName = sampleNameLine.split('sample name =  ',1)[1]
        else:
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            print "echo 'No sample name in %s'" % (errFile.split('.',1)[0] + '.stdout')
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            continue

        #print sampleName
        fileNameLine = sampleLines[3+index]
        fileName = fileNameLine.split('file:',1)[1]
        #print fileName
        jobNumber_1 = errFile.split('ttV_Jan23_',1)[1]
        jobNumber = jobNumber_1.split('.std',1)[0]
        #print jobNumber

        rm_command = "rm %s.std*" % ( errFile.split('.',1)[0] )
        #print rm_command
        print sampleLines[0]
        if "ttV_Jan23" in errFile:
            command = "ttV ssCondor.py %s ttV_Jan23 NA %s 1 %s !>& ttV_%s_%s.log &" % ( sampleName, jobNumber, fileName, sampleName, jobNumber )
            #print command
        

    for outFile in os.popen("find batch_trees/condor_logs/*ttV_Jan23*.stdout -print"):
        sampleLines = open(outFile.strip('\n')).read().splitlines()
        if len(sampleLines) < 3:
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            print "echo 'No sample name in %s'" % outFile
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            continue
        if 'Passed cuts' in sampleLines[len(sampleLines)-2]:
            #print '%s exists, but .stdout file has "Passed cuts"' % errFile 
            continue

        index = 0
        if '.crc.nd.edu' in sampleLines[0]:
            index = 1
        sampleNameLine = sampleLines[2+index]
        if 'sample name =' in sampleNameLine:
            sampleName = sampleNameLine.split('sample name =  ',1)[1]
        else:
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            print "echo 'No sample name in %s'" % outFile
            print "echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'"
            continue

        #print sampleName
        fileNameLine = sampleLines[3+index]
        fileName = fileNameLine.split('file:',1)[1]
        #print fileName
        jobNumber_1 = outFile.split('ttV_Jan23_',1)[1]
        jobNumber = jobNumber_1.split('.std',1)[0]
        #print jobNumber

        rm_command = "rm %s.std*" % ( outFile.split('.',1)[0] )
        print rm_command
        #print sampleLines[0]
        if "ttV_Jan23" in outFile:
            command = "ttV ssCondor.py %s ttV_Jan23 NA %s 1 %s !>& ttV_%s_%s.log &" % ( sampleName, jobNumber, fileName, sampleName, jobNumber )
            print command

        
if __name__ == '__main__':
        main()
        

