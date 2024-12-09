class PAF():
    '''
    Class for dealing with files in PAF format. 
    '''

    def __init__(self):
        '''
        Initialize empty class instance
        '''
        self.alignments = []

    def read(self, filePath):
        '''
        PAF file reader. Read and store data line objects into a list:
        '''
        pafFile = open(filePath)

        # For line in the file
        for line in pafFile:
            
            # Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            fields = line.split() 

            ## Skip truncated lines
            if len(fields) < 12:
                print('[ERROR] Truncated line: ' + str(len(fields)) + ' fields')
                continue

            line = PAF_alignment(fields)
            self.alignments.append(line)

    def sortByLen(self):
        '''
        Sort alignments by query alignment length in descending order
        '''
        sortedAlignments = sorted(self.alignments, key=lambda x: x.alignment_len(), reverse=True)
        return sortedAlignments

    ## [SR CHANGES]
    def sortNbMatches(self):
        '''
        Sort alignments by query alignment length
        '''
        sortedAlignments = sorted(self.alignments, key=lambda x: x.nbMatches, reverse=True)
        return sortedAlignments

    def chain(self, maxDist, maxPercOverlap):
        '''
        Chain PAF alignments based on alignment complementariety

        Input:
            1. maxDist: maximum distance between both ranges 
            2. maxPercOverlap: maximum percentage of overlap between ranges

        Output:
            1. chain: PAF_chain object instance
        '''
        ## 1. Sort alignments by decreasing query alignment length 
        sortedAlignments = self.sortByLen()

        ## 2. Pick longest alignment and initiate chain
        longest = sortedAlignments[0]
        chain = PAF_chain([longest])
        
        # remove alignment 
        del sortedAlignments[0]

        roundCounter = 1

        ## 3. Attemp to extend the chain with complementary alignments
        while True:

            # START ALIGNMENT CHAIN EXTENSION ROUND
            # Initialize boolean as not complementary alignment found
            complBool = False

            # Go through all the available alignments 
            for index, alignment in enumerate(sortedAlignments):

                ## Assess if alignment complementary to the chain
                chainBeg, chainEnd = chain.interval()
                complBool, orientation = gRanges.complementary(chainBeg, chainEnd, alignment.qBeg, alignment.qEnd, maxDist, maxPercOverlap)

                ## Complementary alignment found 
                if complBool:

                    ## Add alignment to the chain 
                    # a) Add to the chain begin
                    if orientation == "LEFT":
                        chain.alignments.insert(0,alignment)

                    # b) Add to the chain end
                    else:
                        chain.alignments.append(alignment)

                    ## Remove from list
                    del sortedAlignments[index]

                    ## Stop once complementary found
                    break

            roundCounter += 1

            # STOP CHAIN EXTENSION IF: 
            # a) No complementary alignment found in the last round OR
            # b) Mo alignments left
            if complBool == False or not sortedAlignments:
                break

        return chain

    def hits2dict(self):
        '''
        Reorganize hits into a dictionary
        '''

        hitsDict = {}

        # For each hit
        for hit in self.alignments:

            # Initialize list
            if hit.qName not in hitsDict:
                hitsDict[hit.qName] = []
            
            # Add hit to list
            hitsDict[hit.qName].append(hit)
    
        return hitsDict

class PAF_alignment():
    '''
    PAF entry class 
    '''
    number = 0 # Number of instances

    def __init__(self, fields):
        '''
        Initialize paf line
        '''
        PAF_alignment.number += 1 # Update instances counter
        self.id = 'PAF_alignment_' + str(PAF_alignment.number)
        self.qName = str(fields[0])
        self.qLen = int(fields[1])
        self.qBeg = int(fields[2])
        self.qEnd = int(fields[3])
        self.strand = str(fields[4])
        self.tName = str(fields[5])
        self.tLen = int(fields[6])
        self.tBeg = int(fields[7])
        self.tEnd = int(fields[8])
        self.nbMatches = int(fields[9])
        self.blockLen = int(fields[10])
        self.MAPQ = int(fields[11])   

    def alignmentLen(self):
        '''
        Compute the query alignment length
        '''

        return self.qEnd - self.qBeg

    def alignmentPerc(self):
        '''
        Compute the query alignment length percentage
        '''
        percLen = float(self.alignmentLen()) / self.qLen * 100

        return percLen

class PAF_chain():
    '''    
    Chain of complementary PAF alignments  
    '''

    def __init__(self, alignments):
        '''
        Initialize chain instance. 
        
        Input:
            1. alignments. List of PAF_alignment instances
        '''
        self.alignments = alignments

    def interval(self):
        '''
        Return query interval covered by the chain
        '''
        firstAlignment = self.alignments[0]
        lastAlignment = self.alignments[-1]
        
        return firstAlignment.qBeg, lastAlignment.qEnd

    def interval_template(self):
        '''
        Return query interval covered by the chain
        '''
        firstAlignment = self.alignments[0]
        lastAlignment = self.alignments[-1]
        
        return firstAlignment.tBeg, lastAlignment.tEnd

    def perc_query_covered(self):
        '''
        Compute the percentage of the query sequence covered by the chain of alignments
        '''
        # a) No alignments available
        if len(self.alignments) == 0:
            percCovered = 0

        # b) Alignments available
        else:

            ##Compute the number of bases covered
            beg, end = self.interval()
            alignmentLen = end - beg

            ## Compute the percentage of bases covered
            percCovered = float(alignmentLen)/self.alignments[0].qLen*100

        return percCovered

    def nb_templates(self):
        '''
        Compute the number of different templates in the chain
        '''
        templates = list(set([alignment.tName for alignment in self.alignments]))
        nbTemplates = len(templates)

        return nbTemplates, templates
