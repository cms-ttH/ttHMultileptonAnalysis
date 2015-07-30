import os

list_name = 'ttV_lists.txt'
label = 'ttV_v1'

def main ():
    listOfLists = open(list_name).read().splitlines()[2:]

    for line in listOfLists:
        try:
            #print line
            list_files = open('../../listsForSkims2012_53x_v3_hadoop/%s.list' % (line)).read().splitlines()[0:]
            input_lists = len(list_files)
            output_files =  len([name for name in os.listdir('batch_trees/%s_lepEff_March18_lepVars' % (line)) if 'root' in name])
            if input_lists != output_files:
                print '%d input lists and %d output files in %s' % (input_lists, output_files, line)
            #print input_lists
            #print output_files
            #else:
                #print '%s is fine' % (line)
        except:
            print 'Sample %s not found' % (line)
        
#         command1 = 'wc -l ../../listsForSkims2012_53x_v3_hadoop/%s.list' % (line)
#         command2 = 'ls batch_trees/%s_lepMVAStudies_v0/*.root | wc -l' % (line)
#         print command1
#         print command2
#         #exec command
            
if __name__ == '__main__':
        main()


        

