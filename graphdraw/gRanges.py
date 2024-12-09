'''
Module 'gRanges' - Contains functions to do operations with genomic ranger or coordinates (i.e. overlaps...)
'''

def overlap_old(start_a, end_a, start_b, end_b):
	"""
	This returns whether there is an overlap between the two or not
	return a Boolean and the overlap length
	"""
	if (start_a < start_b) and (start_b < end_a < end_b):  # overlap from the left
		return True, (end_a - start_b)
	if (end_a > end_b) and (start_b < start_a < end_b):  # overlap from right
		return True, (end_b - start_a)
	else:
		return False, 0


def overlap(begA, endA, begB, endB):
	"""
	check if two ranges overlap and return the overlap length
	"""
	maxBeg = max([begA, begB])
	minEnd = min([endA, endB])

	if maxBeg <= minEnd:
		boolean = True
		overlaplen = minEnd - maxBeg + 1
	else:
		boolean = False
		overlaplen = 0
	return boolean, overlaplen


def complementary(begA, endA, begB, endB, maxDist, maxPercOverlap):
	'''
	Check if two intervals are complementary. Two ranges are considered complementary when:
	1) Ranges located less than X bp of distance AND
	2) Overlapping region does not span X% or more of any of the intervals

	Input:
		1. begA: begin position of range A
		2. endA: end position of range A
		3. begB: begin position of range B
		4. endB: end position of range B
		5. maxDist: maximum distance between both ranges 
		6. maxPercOverlap: maximum percentage of overlap between ranges

	Output:
		1. boolean: True (intervals are complementary) or False (intervals not complementary)
		2. orientation: 'LEFT' (B on the left of A), 'RIGHT' (B on the right of A) or None (not complementary)
	'''  
	## 1. Redefine intervals by adding the maximum distance to their begin and end coordinates
	newBegA = int(begA) - maxDist
	newEndA = int(endA) + maxDist
	newBegB = int(begB) - maxDist
	newEndB = int(endB) + maxDist

	## 2. Assess if redefined intervals do overlap
	boolean = overlap(newBegA, newEndA, newBegB, newEndB)[0]
	
	#A) No overlap
	if not boolean:
		orientation = None
	
	# B) Overlap
	else:

		# 3. Discard those cases with an overlapping region spanning X% or more of at least one of the intervals
		# --------A---------
		#			 <---> Overlapping region
		#			  -------B-------
		overlapLen = overlap(begA, endA, begB, endB)[1]   #Compute real degree of overlap

		#For each interval, compute its % overlapping the other interval
		lenA = (endA - begA) + 1
		lenB = (endB - begB) + 1
		percA = (overlapLen / lenA) * 100
		percB = (overlapLen / lenB) * 100

		# a) One of the intervals overlaps the other by more than the maximum allowed % cutoff 
		if (percA > maxPercOverlap) or (percB > maxPercOverlap):
			boolean = False
			orientation = None

		#b) Overlapping region does not span X% or more of any of the intervals
		else:
			boolean = True

			##Determine complementariety orientation (use A as reference)
			#a) B located on the left of A
			#				  begA <------A------> 
			# begB <------B------>
			if begA > begB:
				orientation = 'LEFT'

			# b) B located on the right of A
			# begA <------A------>
			#				 begB <------B------> 
			else:
				orientation = 'RIGHT'

	return boolean, orientation
