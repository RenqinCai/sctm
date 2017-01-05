class _ChildDoc:
	def __init__(self):
		self.m_ID = -1
		self.m_wordNum = 0
		self.m_inferWordList = []
		self.m_perplexityWordList = []

class _ParentDoc:
	def __init__(self):
		self.m_ID = -1
		self.m_childList = []

class _Corpus:
	def __init__(self):
		self.m_parentList = []
		self.m_parentMap = {}

def readComment(file, corpusObj, proportionThreshold):
	f = open(file)

	rawLine = f.readline()
	line = rawLine.strip().split("\t")
	parentNum = int(line[0])

	parentIndex = 0

	while True:
		rawLine = f.readline()

		if not rawLine:
			break
		line = rawLine.strip().split("\t")

		parentIndex += 1

		if parentIndex > parentNum:
			break

		parentObj = _ParentDoc()
		parentObj.m_ID = parentIndex

		corpusObj.m_parentMap.setdefault(parentIndex, parentObj)
		corpusObj.m_parentList.append(parentIndex)

		commentNum = int(line[0])
		for i in range(commentNum):
			rawLine = f.readline()
			line = rawLine.strip().split("\t")

			childObj = _ChildDoc()
			childObj.m_ID = i
			parentObj.m_childList.append(childObj)

			lineLen = len(line)

			wordNum = int(line[0])
			inferWordLen = int(wordNum*proportionThreshold)

			childObj.m_wordNum = wordNum

			for wordIndex in range(1, inferWordLen+1):
				word = line[wordIndex]
				childObj.m_inferWordList.append(word)

			# if inferWordLen == wordNum:
			# 	print "no word for perplexity"

			if inferWordLen < 1:
				print "no word for infer"

			for wordIndex in range(inferWordLen+1, lineLen):
				word = line[wordIndex]
				childObj.m_perplexityWordList.append(word)

	f.close()

def outputComment(file1, file2, corpusObj):
	f1 = open(file1, "w")
	f2 = open(file2, "w")

	parentNum = len(corpusObj.m_parentList)
	f1.write(str(parentNum)+"\n")

	f2.write(str(parentNum)+"\n")

	for parentID in corpusObj.m_parentList:
		parentObj = corpusObj.m_parentMap[parentID] 

		childNum = len(parentObj.m_childList)
		f1.write(str(childNum)+"\n")
		f2.write(str(childNum)+"\n")

		for childObj in parentObj.m_childList:
			inferLen = len(childObj.m_inferWordList)
			f1.write(str(inferLen)+"\t")

			for word in childObj.m_inferWordList:
				f1.write(word+"\t")

			f1.write("\n")

			perplexityLen = len(childObj.m_perplexityWordList)
			f2.write(str(perplexityLen)+"\t")

			for word in childObj.m_perplexityWordList:
				f2.write(word+"\t")

			f2.write("\n")


	f1.close()
	f2.close()


cv = 10
for i in range(cv):
	corpusObj = _Corpus()
	proportionThreshold = 0
	inputFile = "../input/yahoo_100/testFolder"+str(i)+"/cbagf.AT.txt"
	outputInferFile = "../input/yahoo_100/testFolder"+str(i)+"/cbagf_infer.AT.txt"
	outputPerplexityFile = "../input/yahoo_100/testFolder"+str(i)+"/cbagf_perplexity.AT.txt"

	readComment(inputFile, corpusObj, proportionThreshold)

	outputComment(outputInferFile, outputPerplexityFile, corpusObj)
