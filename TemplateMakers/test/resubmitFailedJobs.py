import os

def main ():

    for errFile in os.popen("find batch_trees/condor_logs/*lepEff_March20_lepVars*.stderr -size +0k -print"):
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
        jobNumber_1 = errFile.split('lepEff_March20_lepVars_',1)[1]
        jobNumber = jobNumber_1.split('.std',1)[0]
        #print jobNumber

        rm_command = "rm %s.std*" % ( errFile.split('.',1)[0] )
        print rm_command
        print 'Failed on node %s' % sampleLines[0]
        if "lepEff_March20_lepVars" in errFile:
            command = "lepEff ssCondor.py %s lepEff_March20_lepVars NA %s 1 %s !>& lepEff_%s_%s.log &" % ( sampleName, jobNumber, fileName, sampleName, jobNumber )
            print command
        

    for outFile in os.popen("find batch_trees/condor_logs/*lepEff_March20_lepVars*.stdout -print"):
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
        jobNumber_1 = outFile.split('lepEff_March20_lepVars_',1)[1]
        jobNumber = jobNumber_1.split('.std',1)[0]
        #print jobNumber

        rm_command = "rm %s.std*" % ( outFile.split('.',1)[0] )
        print rm_command
        print 'Failed on node %s' % sampleLines[0]
        if "lepEff_March20_lepVars" in outFile:
            command = "lepEff ssCondor.py %s lepEff_March20_lepVars NA %s 1 %s !>& lepEff_%s_%s.log &" % ( sampleName, jobNumber, fileName, sampleName, jobNumber )
            print command

        
if __name__ == '__main__':
        main()
        

