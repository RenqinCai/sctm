import numpy as np

class _ParentDoc:
	def __init__(self):
		self.m_ID = -1

class _ChildDoc:
	def __init__(self):
		self.m_ID = -1
		self.m_topicProportion = []
		self.m_wordList = []

class _Corpus:
	def __init__(self):
		self.m_topicNum = 0
		self.m_beta = {}
		###childName:childObj
		self.m_childList = []

def readCommentProportion(file, corpusObj):
	f = open(file)

	while True:
		rawLine = f.readline()

		if not rawLine:
			break

		line = rawLine.strip().split(" ")

		parentID = line[0]
		childNum = int(line[1])

		for i in range(childNum):
			cDoc = _ChildDoc()
			cDoc.m_ID = i

			childName = parentID+"_"+str(i)

			rawLine = f.readline()
			line = rawLine.strip().split(" ")

			lineLen = len(line)
			for i in range(lineLen):
				topicPro = float(line[i])
				cDoc.m_topicProportion.append(topicPro)

			corpusObj.m_childList.append(cDoc)

		rawLine = f.readline()

	f.close()

def readBeta(file, corpusObj):
	f = open(file)

	rawLine = f.readline()

	line = rawLine.strip().split(" ")

	topicNum = int(line[0])
	wordNum = int(line[1])

	corpusObj.m_topicNum = topicNum

	for i in range(topicNum):
		rawLine = f.readline()
		line = rawLine.strip().split(" ")

		lineLen = len(line)

		corpusObj.m_beta.setdefault(i, [])
		topicWordProbList = []

		for wordIndex in range(lineLen):
			topicWordProb = float(line[wordIndex])
			topicWordProbList.append(topicWordProb)

		corpusObj.m_beta[i] = topicWordProbList

	f.close()

def readCommentWord(file, corpusObj):
	f = open(file)

	rawLine = f.readline()
	line = rawLine.strip().split("\t")

	parentNum = int(line[0])

	parentIndex = 0

	commentIndex = 0

	while True:
		rawLine = f.readline()

		if not rawLine:
			break

		line = rawLine.strip().split("\t")


		parentIndex += 1
		if parentIndex > parentNum:
			break

		commentNum = int(line[0])
		for i in range(commentNum):
			rawLine = f.readline()
			line = rawLine.strip().split("\t")
			# print line
			childObj = corpusObj.m_childList[commentIndex]

			lineLen = len(line)
			wordNum = int(line[0])

			for wordIndex in range(1, lineLen):
				word = int(line[wordIndex])
				childObj.m_wordList.append(word)

			commentIndex += 1

	print commentIndex

	f.close()

def computePerplexity(corpusObj):
	beta = {}
	beta = corpusObj.m_beta
	topicNum = corpusObj.m_topicNum

	totalLikelihood = 0
	totalWordNum = 0
	totalChildNum = 0

	for childObj in corpusObj.m_childList:
		wordList = childObj.m_wordList
		totalWordNum += len(wordList)
		if len(wordList) == 0:
			print "no word in comment"
		totalChildNum += 1

		for word in wordList:
			likelihood = 0
			for topicIndex in range(topicNum):
				theta = childObj.m_topicProportion[topicIndex]
				betaProb = beta[topicIndex][word]

				likelihood += theta*betaProb

			# print likelihood
			if likelihood < 1e-20:
				print "small value likelihood"
				likelihood += 1e-20
			totalLikelihood += np.log(likelihood)

	print("likelihood\t"+str(totalLikelihood))
	print("wordNum \t"+str(totalWordNum))
	print("childNum\t"+str(totalChildNum))
	perplexity = np.exp(-totalLikelihood*1.0/totalWordNum)
	print("perplexity\t"+str(perplexity))

commentProportionFile = "../output/tech_100/0/topic15/y_dist_test.txt"
betaFile = "../output/tech_100/0/topic15/beta"

commentWordFile = "../input/tech_100/testFolder0/cbagf_perplexity.AT.txt"

corpusObj = _Corpus()

readCommentProportion(commentProportionFile, corpusObj)
readBeta(betaFile, corpusObj)
readCommentWord(commentWordFile, corpusObj)
computePerplexity(corpusObj)


