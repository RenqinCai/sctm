import numpy as np

class _parentDoc:
	def __init__(self, ID):
		self.m_ID = ID
		###sentenceID:sentenceObj
		self.m_sentenceList = []
		self.m_childList = []

class _sentence:
	def __init__(self, ID):
		self.m_ID = ID
		self.m_wordsList = []

		###childID:likelihoodValue
		self.m_childLikelihoodMap = {}

class _childDoc:
	def __init__(self, ID):
		self.m_ID = ID
		self.m_parentDocID = ""
		self.m_topicProportion = []

class _corpus:
	def __init__(self):
		self.m_parentMap = {}
		self.m_childMap = {}
		self.m_beta = {}
		self.m_topicNum = 0
		self.m_wordNum = 0

def readSelectedStns(selectedStnFile, stnWordFile, corpusObj):
	selectedStnF = open(selectedStnFile)

	for rawLine in selectedStnF:
		line = rawLine.strip().split("\t")

		parentID = line[0]
		stnNum = int(line[1])

		pDocObj = _parentDoc(parentID)
		corpusObj.m_parentMap.setdefault(parentID, pDocObj)

		lineLen = len(line)

		for lineIndex in range(2, lineLen):
			stnIndex = line[lineIndex]

			stnObj = _sentence(stnIndex)
			pDocObj.m_sentenceList.append(stnObj)


	selectedStnF.close()

	stnWordF = open(stnWordFile)

	rawLine = stnWordF.readline()
	line = rawLine.strip().split("\t")

	parentNum = int(line[0])

	parentIndex = 0
	stnIndex = 0

	while True:
		rawLine = stnWordF.readline()

		if not rawLine:
			break

		line = rawLine.strip().split("\t")

		parentObj = corpusObj.m_parentMap[str(parentIndex+1)]
		parentIndex += 1

		if parentIndex > parentNum:
			break

		stnNum = int(line[0])

		for stnIndex in range(stnNum):
			rawLine = stnWordF.readline()

			line = rawLine.strip().split("\t")

			stnObj = parentObj.m_sentenceList[stnIndex]

			lineLen = len(line)
			wordNum = int(line[0])

			for wordIndex in range(1, lineLen):
				word = int(line[wordIndex])
				stnObj.m_wordsList.append(word)

	stnWordF.close()

def readSelectedComments(commentFile, corpusObj):
	f = open(commentFile)

	for rawLine in f:
		line = rawLine.strip().split("\t")

		parentID = line[0]
		parentObj = corpusObj.m_parentMap[parentID]

		lineLen = len(line)
		for lineIndex in range(1, lineLen):
			childName = line[lineIndex]
			parentObj.m_childList.append(childName)

	f.close()

def readCommentProportion(proportionFile, corpusObj):
	f = open(proportionFile)

	while True:
		rawLine = f.readline()

		if not rawLine:
			break

		line = rawLine.strip().split(" ")

		parentID = line[0]
		childNum = int(line[1])

		parentObj = corpusObj.m_parentMap[parentID]

		for childIndex in range(childNum):
			cDoc = _childDoc(childIndex)

			childName = parentObj.m_childList[childIndex]

			rawLine = f.readline()
			line = rawLine.strip().split(" ")

			lineLen = len(line)
			for lineIndex in range(lineLen):
				topicPro = float(line[lineIndex])
				cDoc.m_topicProportion.append(topicPro)

			corpusObj.m_childMap.setdefault(childName, cDoc)
			# parentObj.m_childList.append(childName)

		rawLine = f.readline()

	f.close()

def readBeta(file, corpusObj):
	f = open(file)

	rawLine = f.readline()

	line = rawLine.strip().split(" ")

	topicNum = int(line[0])
	wordNum = int(line[1])

	corpusObj.m_topicNum = topicNum

	for topicIndex in range(topicNum):
		rawLine = f.readline()
		line = rawLine.strip().split(" ")

		lineLen = len(line)
		corpusObj.m_beta.setdefault(topicIndex, [])
		topicWordProbList = []

		for wordIndex in range(lineLen):
			topicWordProb = float(line[wordIndex])
			topicWordProbList.append(topicWordProb)

		corpusObj.m_beta[topicIndex] = topicWordProbList

	f.close()

def computeSentenceLikelihood(corpusObj):
	beta = {}
	beta = corpusObj.m_beta

	topicNum = corpusObj.m_topicNum

	for parentID in corpusObj.m_parentMap.keys():
		parentObj = corpusObj.m_parentMap[parentID]

		for stnObj in parentObj.m_sentenceList:
			wordList = stnObj.m_wordsList

			for childName in parentObj.m_childList:
				cDoc = corpusObj.m_childMap[childName]

				childLikelihood = 0

				for word in wordList:
					likelihood = 0

					for topicIndex in range(topicNum):
						theta = cDoc.m_topicProportion[topicIndex]

						betaProb = beta[topicIndex][word]

						likelihood += theta*betaProb

					if likelihood < 1e-20:
						likelihood += 1e-20

					childLikelihood += np.log(likelihood)
				
				stnObj.m_childLikelihoodMap.setdefault(childName, childLikelihood)

def outputTopChild4Stn(topChild4StnFile, corpusObj):

	f = open(topChild4StnFile, "w")

	for parentID in corpusObj.m_parentMap.keys():
		parentObj = corpusObj.m_parentMap[parentID]

		f.write(parentID+"\t")

		stnNum = len(parentObj.m_sentenceList)

		f.write(str(stnNum)+"\n")

		for stnObj in parentObj.m_sentenceList:
			f.write(stnObj.m_ID+"\t")

			for childName in parentObj.m_childList:
				childLikelihood = stnObj.m_childLikelihoodMap[childName]
				f.write(childName+":"+str(childLikelihood)+"\t")

			f.write("\n")

	f.close()


corpusObj = _corpus()

selectedStnFile = "./input/tech_corpus/selected_Stn.txt"
stnWordFile = "./input/tech_corpus/abagf.AT.txt"
selectedCommentFile = "./input/tech_corpus/selected_Comments.txt"

commentProportionFile = "./output/tech_corpus/y_dist.txt"
betaFile = "./output/tech_corpus/beta"

outputFile = "./output/tech_corpus/topChild4Stn_SCTM.txt"

readSelectedStns(selectedStnFile, stnWordFile, corpusObj)
readSelectedComments(selectedCommentFile, corpusObj)

readCommentProportion(commentProportionFile, corpusObj)
readBeta(betaFile, corpusObj)

computeSentenceLikelihood(corpusObj)
outputTopChild4Stn(outputFile, corpusObj)





