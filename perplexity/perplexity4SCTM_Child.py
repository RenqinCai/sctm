import numpy as np

class _ParentDoc:
	def __init__(self):
		self.m_ID = -1
		self.m_wordList = []
		self.m_topicProportion = []

class _ChildDoc:
	def __init__(self):
		self.m_ID = -1
		self.m_topicProportion = []
		self.m_wordList = []
		self.m_parentIndex = 0

class _Corpus:
	def __init__(self):
		self.m_topicNum = 0
		self.m_beta = {}
		###childName:childObj
		self.m_childList = []
		self.m_parentList = []

def readArticleProportion(file, corpusObj):
	f = open(file)

	while True:
		rawLine = f.readline()

		if not rawLine:
			break


		line = rawLine.strip().split(" ")

		parentID = line[0]

		proportionLine = f.readline()
		pDoc = _ParentDoc()
		pDoc.m_ID = parentID

		proportionLine = proportionLine.strip().split(" ")
		topicLen = len(proportionLine)
		for i in range(topicLen):	

			topicProportion = float(proportionLine[i])
			# print "topic Proportion\t",topicProportion
			pDoc.m_topicProportion.append(topicProportion)

		corpusObj.m_parentList.append(pDoc)

	f.close()

def readCommentProportion(file, corpusObj):
	f = open(file)

	while True:
		rawLine = f.readline()

		if not rawLine:
			break

		line = rawLine.strip().split(" ")

		parentID = int(line[0])-1
		childNum = int(line[1])

		for i in range(childNum):
			cDoc = _ChildDoc()
			cDoc.m_ID = i

			# childName = parentID+"_"+str(i)
			cDoc.m_parentIndex = parentID

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


		# parentIndex += 1
		# if parentIndex > parentNum:
		# 	break


		commentNum = int(line[0])
		for i in range(commentNum):
			rawLine = f.readline()
			line = rawLine.strip().split("\t")
			# print line
			childObj = corpusObj.m_childList[commentIndex]

			if childObj.m_parentIndex != parentIndex:
				print "error"

			lineLen = len(line)
			wordNum = int(line[0])

			for wordIndex in range(1, lineLen):
				word = int(line[wordIndex])
				childObj.m_wordList.append(word)

			commentIndex += 1

		parentIndex += 1

	print commentIndex

	f.close()

def readArticleWord(file, corpusObj):
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

		parentObj = corpusObj.m_parentList[parentIndex]
	
		sentenceNum = int(line[0])
		for i in range(sentenceNum):
			rawLine = f.readline()
			line = rawLine.strip().split("\t")
			# print line

			lineLen = len(line)
			wordNum = int(line[0])

			for wordIndex in range(1, lineLen):
				word = int(line[wordIndex])
				parentObj.m_wordList.append(word)

		parentIndex += 1
		if parentIndex > parentNum:
			break
	f.close()

def computePerplexity(corpusObj):
	beta = {}
	beta = corpusObj.m_beta
	topicNum = corpusObj.m_topicNum

	totalLikelihood = 0
	totalWordNum = 0
	totalChildNum = 0

	for childObj in corpusObj.m_childList:
		parentIndex = childObj.m_parentIndex
		parentObj = corpusObj.m_parentList[parentIndex]

		wordList = childObj.m_wordList
		totalWordNum += len(wordList)
		if len(wordList) == 0:
			print "no word in comment"
		totalChildNum += 1

		for word in wordList:
			likelihood = 0
			for topicIndex in range(topicNum-1):
				theta = parentObj.m_topicProportion[topicIndex]
				theta = 1
				betaProb = beta[topicIndex][word]

				likelihood += theta*betaProb

			# print likelihood
			if likelihood < 1e-20:
				print "child small value likelihood"
				likelihood += 1e-20
			totalLikelihood += np.log(likelihood)

	# for parentObj in corpusObj.m_parentList:
	# 	wordList = parentObj.m_wordList
	# 	totalWordNum += len(wordList)

	# 	if len(wordList) == 0:
	# 		print "no word in comment"


	# 	topicProportionLen = len(parentObj.m_topicProportion)
	# 	if topicProportionLen != topicNum -1:
	# 		print "error topic proportion len"

	# 	for word in wordList:
	# 		likelihood = 0
	# 		for topicIndex in range(topicNum-1):
	# 			theta = parentObj.m_topicProportion[topicIndex]

	# 			betaProb = beta[topicIndex][word]

	# 			likelihood += theta*betaProb


	# 		if likelihood < 1e-20:
	# 			print likelihood
	# 			print "parent small value likelihood"
	# 			likelihood += 1e-20

	# 		totalLikelihood += np.log(likelihood)


	print("likelihood\t"+str(totalLikelihood))
	print("wordNum \t"+str(totalWordNum))
	print("childNum\t"+str(totalChildNum))
	perplexity = np.exp(-totalLikelihood*1.0/totalWordNum)
	print("perplexity\t"+str(perplexity))

commentProportionFile = "../output/tech_100/0/topic60/y_dist_test.txt"
articleProportionFile = "../output/tech_100/0/topic60/z_distDoc_test.txt"

betaFile = "../output/tech_100/0/topic60/beta"

commentWordFile = "../input/tech_100/testFolder0/cbagf_perplexity.AT.txt"
# articleWordFile = "../input/tech_80/testFolder0/abagf.AT.txt"

corpusObj = _Corpus()

readArticleProportion(articleProportionFile, corpusObj)
readCommentProportion(commentProportionFile, corpusObj)

readBeta(betaFile, corpusObj)
# readArticleWord(articleWordFile, corpusObj)
readCommentWord(commentWordFile, corpusObj)
computePerplexity(corpusObj)


