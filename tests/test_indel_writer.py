import pytest
from phastSim.phastSim import mutation, mType
from copy import deepcopy



class DummyNode:

    def __init__(self, mutationsList):
        self.mutations = mutationsList
        self.name = "dummy_node"

class WriterTest:


    def __init__(self):
        self.ref = "A" * 2000
        self.refList = list(self.ref)


    def getAbsoluteCoords(self, mut, mutDict, insertionDict):
        """
        refPos - the position on the reference that this mut appears on. 
        offset - if the mutation is nested inside another mutation, or not at the start of a series of mutations, 
        then it may have a non-zero offset, referring to how far along the mutation is inside another already existing mutation. 
        """
        if mut.insertionPos:
            
            refPos = insertionDict[mut.insertionPos]
            for offset, coords in enumerate(mutDict[refPos][2]):
                if (mut.insertionPos, mut.genomePos) == coords:
                    return refPos, offset

        else:
            return mut.genomePos, 0

    def wrapWriterWithMutationDict(writingFunction):
        """
        This function wraps a writing function, by keeping track of a mutDict and insertionDict; the writing function is there
        just to write the output to the file in the correct format. 

        I've done this because the same wrapper with a mutDict and insertionDict is used in multiple places. 

        This whole function is overly complicated and potentially could be quite slow. 
        The aim is to 'collapse' all of the (possibly overlapping / nested) indels and substitutions
        into a concise format. 

        Any advice on how to avoid all this complicated work (either by simplifying the code or deciding on
        a different but still unambiguous output format) is appreciated!

        mutDict = {
            position: (original nuc, mutated nucs, list of the coordinates of the new nuc),

            1: (A, -AATTTG, [(0,1), (2,0), (2,1),(5, 0), (5,1), (5,2), (2,2)])

            123: (ATT, "", [])

            456: (A, G, [(0, 456)])
        }

        insertionDict = {
            2: 1
            5: 1
        }

        (2, 1) -> refPos = 1, offset = 2

        (0, 123) -> refPos = 0, offset = 123
        
        mutation(source = "", target = T, insertionPos = 0, genomePos = 123, [index = 6])

        mutation(source = AAAAA, target = "", insertionPos = 0, genomePos = 234)

        """

        def wrappedFunction(self, node, file, mutDict, insertionDict, **kwargs):

            print(f"MutDict BEFORE updating for mutations in {node.name}**************")
            print(mutDict)
            print(f"insertionDict BEFORE updating for mutations in {node.name}**************")
            print(insertionDict)
            print("Mutation events")
            print([str(x) for x in node.mutations])

            # flatten the mutDict + node.mutations into a more concise format, and save indels in insertionDict
            for m in node.mutations:

                # find the 'absolute' coordinate of the mutation i.e. where it lives relative to the 
                # already existing mutations in the mutDict. 
                refPos, offset = self.getAbsoluteCoords(m, mutDict, insertionDict)

                print(f"***************mutation {(str(m))} on branch {node.name} @ coords {refPos}, {offset}**************")
                print({k:v for k,v in mutDict.items() if abs(1268 - k) < 20})
                print(f"**************mutation {str(m)} on branch {node.name}***************")
                print({k:v for k,v in insertionDict.items() if abs(1268 - k) < 10})


                if m.mType == mType.SUB:

                    if refPos in mutDict:
                        target = mutDict[refPos][1]
                        mutDict[refPos][1] = target[:offset] + m.target + target[offset+1:]

                        # if the substitution has exactly reversed a previous substitution then we can delete it
                        if mutDict[refPos][1] == self.ref[refPos]:
                            del mutDict[refPos]

                    else:
                        mutDict[refPos] = [self.ref[refPos], m.target, [(m.insertionPos, m.genomePos)]]

                elif m.mType == mType.INS:

                    # add the insertion to the insertion dictionary
                    insertionDict[m.index] = refPos

                    if refPos in mutDict:
                        # record the indices i.e. the relative genome coordinates of the existing data
                        indices = mutDict[refPos][2]
                        target = mutDict[refPos][1]
                        

                        # the new data is comprised of: original, new insertion, new indices of genome
                        mutDict[refPos] = [mutDict[refPos][0], 
                                        target[:offset+1] + m.target + target[offset+1:],
                                        indices[:offset+1] + [(m.index, i) for i in range(len(m.target))] + indices[offset+1:]]

                    else:
                        # deal with an edge case where we are at the -1 genome position i.e. right at the start of the genome
                        first_symbol = (self.ref[refPos] if refPos != -1 else "")
                        mutDict[refPos] = [first_symbol,
                                        first_symbol + m.target,
                                        [(m.insertionPos, m.genomePos)] + [(m.index, i) for i in range(len(m.target))]]

                # deletion - this is the most difficult case
                else:

                    # need to deal with the deletion one symbol at a time
                    deletedChars = 0
                    m.deletedPositions = []
                    
                    startedInMiddleOfDictionaryItem = False
                    if offset > 0:
                        startedInMiddleOfDictionaryItem = True

                    # are we appending "-" characters to the end of an existing dictionary item or not?
                    appending = False

                    # dictPos - the dictionary key in mutDict that we are reading from / writing to
                    dictPos = refPos

                    while (deletedChars < len(m.source)):

                        # we may be appending characters to the current deletion one at a time
                        if appending:
                            
                            # if the next character is already in the dictionary, we must stop appending 
                            if refPos in mutDict:
                            
                                appending = False
                                
                                # delete a character if it is non-blank
                                if mutDict[dictPos][1][offset] != "-":
                                    mutDict[dictPos][1] = mutDict[dictPos][1][:offset] + "-" + mutDict[dictPos][1][offset+1:]
                                    deletedChars += 1
                                    m.deletedPositions.append(mutDict[dictPos][2][offset])
                            
                            # otherwise we can continue appending to the deletion 
                            else:
                                
                                mutDict[dictPos][0] += self.ref[refPos]
                                mutDict[dictPos][1] += "-"
                                mutDict[dictPos][2] += [(0, refPos)]
                                deletedChars += 1
                                m.deletedPositions.append((0, refPos))
                        
                        # in this case we are working on an existing dictionary item, or have just finished doing so
                        else:
                            
                            # a position not ever seen before, we must add it to the dictionary and 
                            # can append to it
                            if not dictPos in mutDict:
                                
                                mutDict[dictPos] = [self.ref[dictPos], "-", [(0, dictPos)]]
                                deletedChars += 1
                                m.deletedPositions.append((0, dictPos))
                                appending = True
                                
                            # delete a character if it is non-blank
                            else: 
                                if mutDict[dictPos][1][offset] != "-":
                                    mutDict[dictPos][1] = mutDict[dictPos][1][:offset] + "-" + mutDict[dictPos][1][offset+1:]
                                    deletedChars += 1
                                    m.deletedPositions.append(mutDict[dictPos][2][offset])
                        
                        
                        # now move to the next character
                        # increment the reference position if we pass it
                        if (mutDict[dictPos][2][offset] == (0, refPos)):
                            refPos += 1  
                            
                        # we will change to a new dictionary item if we are not in append mode or are not going to be
                        offset += 1
                        if offset == len(mutDict[dictPos][2]):
                            if (not appending) or (refPos in mutDict):
                                
                                if startedInMiddleOfDictionaryItem:
                                    refPos += 1
                                    startedInMiddleOfDictionaryItem = False
                                appending = False
                                dictPos = refPos
                                offset = 0

            print(f"MutDict DURING printing for {node.name}**************")
            print(mutDict)
            print(f"insertionDict DURING printing for {node.name}**************")
            print(insertionDict)

            # call the writing function - which will write to the file or pass data to the child nodes as appropriate
            writingFunction(self, node, file, mutDict, insertionDict, **kwargs)

            # de-update the list so it can be used by siblings etc.
            # we need to do everything in the first part but backwards - terrific
            for m in reversed(node.mutations):

                refPos, offset = self.getAbsoluteCoords(m, mutDict, insertionDict)

                if m.mType == mType.SUB:

                    if refPos in mutDict:
                        target = mutDict[refPos][1]
                        mutDict[refPos][1] = target[:offset] + m.source + target[offset+1:]

                        # if the substitution has exactly reversed a previous substitution then we can delete it
                        if mutDict[refPos][1] == self.ref[refPos]:
                            del mutDict[refPos]

                    else:
                        mutDict[refPos] = [self.ref[refPos], m.source, [(m.insertionPos, m.genomePos)]]

                elif m.mType == mType.INS:

                    # de-update an insertion by removing it from the insertion dictionary 
                    # and also removing any reference to it in the mutDict
                    del insertionDict[m.index]
                    insertionLength = len(m.target)

                    # remove the insertion
                    mutDict[refPos][1] = mutDict[refPos][1][:offset + 1] + mutDict[refPos][1][offset + insertionLength + 1:]
                    mutDict[refPos][2] = mutDict[refPos][2][:offset + 1] + mutDict[refPos][2][offset + insertionLength + 1:]

                    # if the de-updated mutation is now identical to the reference genome then delete the dictionary item
                    if mutDict[refPos][0] == mutDict[refPos][1]:
                        del mutDict[refPos]

                # deletion - this is the more difficult case - we need to put back everything that had been deleted. 
                # this may mean putting back an arbitrary string of insertions and substitutions
                else:


                    # need to deal with the deletion one symbol at a time
                    deletedChars = 0
                    counter = 0

                    while (deletedChars < len(m.source)):
                        # a deletion may have skipped characters (other nested deletions)
                        # we need to make sure we only 'un-delete' the correct characters
                        if mutDict[refPos][2][offset] == m.deletedPositions[deletedChars]:
                            
                            # insert back symbols in the target mutation
                            mutDict[refPos][1] = mutDict[refPos][1][:offset] + m.source[deletedChars] + mutDict[refPos][1][offset+1:]
                            deletedChars += 1

                        # move to next symbol
                        if mutDict[refPos][2][offset] == (0, refPos + counter):
                            counter += 1
                        offset += 1
                        if offset == len(mutDict[refPos][1]):
                            offset = 0
                            
                            # remove 'null' deletions
                            if mutDict[refPos][0] == mutDict[refPos][1]:
                                del mutDict[refPos]

                            refPos += counter
                            counter = 0

            print(f"MutDict AFTER de-updating for mutations in {node.name}**************")
            print(mutDict)
            print(f"insertionDict AFTER de-updating for mutations in {node.name}**************")
            print(insertionDict)

        return wrappedFunction

    @wrapWriterWithMutationDict
    def null_function(*args, **kwargs):
        pass


def apply_general_test(testCase):

    mutDict, insertionDict, mutationList = testCase
    node = DummyNode(mutationsList=mutationList)

    mutDictBefore = deepcopy(mutDict)
    insertionDictBefore = deepcopy(insertionDict)

    writer = WriterTest()

    writer.null_function(mutDict=mutDict, insertionDict=insertionDict, file=None, node=node)

    assert(mutDict == mutDictBefore)
    assert(insertionDict == insertionDictBefore)


def test_substitution_on_genome():
    
    testCase = ({
        1267: ["A", "ACGTGCTAAC", [(0, 1267)] + [(70, x) for x in range(0, 9)]],
        1276: ["A", "AGCACTTTCTGT", [(0, 1276)] + [(54, x) for x in range(0, 11)]]
    },
    {54: 1276, 70: 1267},
    [mutation(mType=mType.SUB, genomePos=1, insertionPos=0, source="A", target="G")])

    apply_general_test(testCase)


def test_substitution_inside_indel():
    
    testCase = ({
        1267: ["A", "ACGTGCTAAC", [(0, 1267)] + [(70, x) for x in range(0, 9)]],
        1276: ["A", "AGCACTTTCTGT", [(0, 1276)] + [(54, x) for x in range(0, 11)]]
    },
    {54: 1276, 70: 1267},
    [mutation(mType=mType.SUB, genomePos=5, insertionPos=70, source="T", target="G")])

    apply_general_test(testCase)


def test_insertion_on_genome():

    testCase = ({
        1267: ["A", "ACGTGCTAAC", [(0, 1267)] + [(70, x) for x in range(0, 9)]],
        1276: ["A", "AGCACTTTCTGT", [(0, 1276)] + [(54, x) for x in range(0, 11)]]
    },
    {54: 1276, 70: 1267},
    [mutation(mType=mType.INS, genomePos=5, insertionPos=0, source="", target="GCGGGCGC", index=1)])

    apply_general_test(testCase)

def test_insertion_inside_insertion():

    testCase = ({
        1267: ["A", "ACGTGCTAAC", [(0, 1267)] + [(70, x) for x in range(0, 9)]],
        1276: ["A", "AGCACTTTCTGT", [(0, 1276)] + [(54, x) for x in range(0, 11)]]
    },
    {54: 1276, 70: 1267},
    [mutation(mType=mType.INS, genomePos=5, insertionPos=70, source="", target="GCGGGCGC", index=1)])

    apply_general_test(testCase)


def test_deletion_on_genome():

    testCase = ({
        1267: ["A", "ACGTGCTAAC", [(0, 1267)] + [(70, x) for x in range(0, 9)]],
        1276: ["A", "AGCACTTTCTGT", [(0, 1276)] + [(54, x) for x in range(0, 11)]]
    },
    {54: 1276, 70: 1267},
    [mutation(mType=mType.DEL, genomePos=5, insertionPos=0, source="AAAAAAAAA", target="")])

    apply_general_test(testCase)


def test_deletion_overlapping_starting_from_indel():

    testCase = ({
        1267: ["A", "ACGTGCTAAC", [(0, 1267)] + [(70, x) for x in range(0, 9)]],
        1276: ["A", "AGCACTTTCTGT", [(0, 1276)] + [(54, x) for x in range(0, 11)]]
    },
    {54: 1276, 70: 1267},
    [mutation(mType=mType.DEL, genomePos=5, insertionPos=70, source="TAACAAAAAAAAAGC", target="")])

    apply_general_test(testCase)

def test_multiple_deletions_together():

    testCase = ({
        1406: ['A', '----', [(0, 1406), (23, 0), (23, 1), (23, 2)]], 
        1405: ['A', '-', [(0, 1405)]], 
        1401: ['AAAA', '----', [(0, 1401), (0, 1402), (0, 1403), (0, 1404)]],
        1408: ['AAA', '---', [(0, 1408), (0, 1409), (0, 1410)]]
    },
    {23: 1406},
    [mutation(mType=mType.DEL, genomePos=1399, insertionPos=0, source="AAAAAAAAAAAAAA", target="")])

    apply_general_test(testCase)

def test_multiple_deletions_together_slightly_different():

    testCase = ({
        1406: ['A', '-CTG', [(0, 1406), (23, 0), (23, 1), (23, 2)]], 
        1401: ['AAAA', '----', [(0, 1401), (0, 1402), (0, 1403), (0, 1404)]],
        1408: ['AAA', '---', [(0, 1408), (0, 1409), (0, 1410)]]
    },
    {23: 1406},
    [mutation(mType=mType.DEL, genomePos=1399, insertionPos=0, source="AAACTGAAAAAAAA", target="")])

    apply_general_test(testCase)
