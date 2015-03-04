/*
 * coveredIntervals.cpp
 *
 *  Created on: 01.10.2014
 *      Author: AlexanderDilthey
 */

#include "coveredIntervals.h"

#include <assert.h>

namespace GraphAlignerUnique {

coveredIntervals::coveredIntervals() {
	// TODO Auto-generated constructor stub
	intervalTree_upToDate = false;
}

bool coveredIntervals::externalIntervalCovered(std::string regionID, int from, int to) const
{
	std::set<oneInterval*> S = externalIntervalCoveredBy(regionID, from, to);
	assert(S.size() <= 1);
	return (S.size() > 0);
}

void coveredIntervals::buildIntervalTrees()
{
	validate();

	for(std::map<std::string, IntervalTree<oneInterval*>*>::iterator treeIt = intervalTrees.begin(); treeIt != intervalTrees.end(); treeIt++)
	{
		IntervalTree<oneInterval*>* tree = treeIt->second;
		delete(tree);
	}

	std::map<std::string, std::vector<Interval<oneInterval*>>> intervals_by_region;

	for(std::set<oneInterval*>::iterator intervalIt = intervals.begin(); intervalIt != intervals.end(); intervalIt++)
	{
		oneInterval* interval = *intervalIt;

		std::string regionID = interval->regionID;
		int from = interval->from;
		int to = interval->to;

		Interval<oneInterval*> intervalForTree(from, to, interval);

		intervals_by_region[regionID].push_back(intervalForTree);
	}

	for(std::map<std::string, std::vector<Interval<oneInterval*>>>::iterator regionIt = intervals_by_region.begin(); regionIt != intervals_by_region.end(); regionIt++)
	{
		std::string regionID = regionIt->first;
		std::vector<Interval<oneInterval*>> intervals = regionIt->second;

		IntervalTree<oneInterval*>* tree = new IntervalTree<oneInterval*>(intervals);
		intervalTrees[regionID] = tree;
	}

	intervalTree_upToDate = true;


}

bool coveredIntervals::isThere_overlap_with_givenInterval(std::string regionID, int from, int to)
{
	std::set<oneInterval*> overlap = intervals_overlapping_with_givenInterval(regionID, from, to);
	return (overlap.size() > 0);
}

std::set<oneInterval*> coveredIntervals::intervals_overlapping_with_givenInterval(std::string regionID, int from, int to)
{
	if(! intervalTree_upToDate)
	{
		buildIntervalTrees();
	}

	std::set<oneInterval*> forReturn;
	std::vector<Interval<oneInterval*> > R;
	if(intervalTrees.count(regionID))
	{
		intervalTrees.at(regionID)->findOverlapping(from, to, R);
		for(unsigned int rI = 0; rI < R.size(); rI++)
		{
			forReturn.insert(R.at(rI).value);
		}
	}

	return forReturn;
}

std::set<oneInterval*> coveredIntervals::externalIntervalCoveredBy(std::string regionID, int from, int to) const
{
	// implemented in this way because this might be faster if properly done (i.e.intervals_covering_point fast)
	std::set<oneInterval*> candidates1 = intervals_covering_point(regionID, from);
	std::set<oneInterval*> candidates2 = intervals_covering_point(regionID, to);

	std::set<oneInterval*> candidates;
	candidates.insert(candidates1.begin(), candidates1.end());
	candidates.insert(candidates2.begin(), candidates2.end());

	std::set<oneInterval*> forReturn;
	for(std::set<oneInterval*>::iterator iIt = candidates.begin(); iIt != candidates.end(); iIt++)
	{
		if(candidates1.count(*iIt) && candidates2.count(*iIt))
		{
			oneInterval* coveringInterval = *iIt;
			assert((coveringInterval->from <= from) && (coveringInterval->to >= from));
			assert((coveringInterval->from <= to) && (coveringInterval->to >= to));

			forReturn.insert(coveringInterval);
		}
	}

	return forReturn;
}

void coveredIntervals::addPoint(std::string regionID, int position)
{
	intervalTree_upToDate = false;

	if(intervals_covering_point(regionID, position).size() == 0)
	{
		int p_minus1 = position - 1;
		int p_plus1 = position + 1;

		std::set<oneInterval*> right_extensible_intervals = intervals_covering_point(regionID, p_minus1);
		std::set<oneInterval*> left_extensible_intervals = intervals_covering_point(regionID, p_plus1);

		assert(right_extensible_intervals.size() <= 1);
		assert(left_extensible_intervals.size() <= 1);

		if((right_extensible_intervals.size() == 1) && (left_extensible_intervals.size() == 1))
		{
			oneInterval* leftInterval = *(right_extensible_intervals.begin());
			oneInterval* rightInterval = *(left_extensible_intervals.begin());

			assert(leftInterval->to < rightInterval->from);

			leftInterval->to = rightInterval->to;

			removeIntervalFromIndices(rightInterval);
			intervals.erase(rightInterval);

			delete(rightInterval);

		}
		else if(right_extensible_intervals.size() == 1)
		{
			assert(left_extensible_intervals.size() == 0);

			oneInterval* toExtend = *(right_extensible_intervals.begin());
			assert(toExtend->to == (position - 1));
			removeIntervalFromIndices(toExtend);
			toExtend->to = position;
			addIntervalToIndices(toExtend);
		}
		else if(left_extensible_intervals.size() == 1)
		{
			assert(right_extensible_intervals.size() == 0);

			oneInterval* toExtend = *(left_extensible_intervals.begin());
			assert(toExtend->from == (position + 1));
			removeIntervalFromIndices(toExtend);
			toExtend->from = position;
			addIntervalToIndices(toExtend);
		}
		else
		{
			assert((right_extensible_intervals.size() == 0) && (left_extensible_intervals.size() == 0));

			oneInterval* newInterval = new oneInterval;
			newInterval->regionID = regionID;
			newInterval->from = position;
			newInterval->to = position;

			intervals.insert(newInterval);
			addIntervalToIndices(newInterval);
		}
	}
}

void coveredIntervals::removeIntervalFromIndices(oneInterval* i)
{
	intervalTree_upToDate = false;

	assert(intervals_start.at(i->regionID).count(i->from));
	intervals_start.at(i->regionID).erase(i->from);
	assert(intervals_stop.at(i->regionID).count(i->to));
	intervals_stop.at(i->regionID).erase(i->to);
}

void coveredIntervals::addIntervalToIndices(oneInterval* i)
{
	intervalTree_upToDate = false;

	std::string regionID = i->regionID;
	assert((intervals_start.count(regionID) == 0) || (intervals_start.at(regionID).count(i->from) == 0));
	intervals_start[i->regionID][i->from] = i;
	assert((intervals_stop.count(i->regionID) == 0) || (intervals_stop.at(i->regionID).count(i->to) == 0));
	intervals_stop[i->regionID][i->to] = i;
}    

size_t coveredIntervals::getNumIntervals() const
{
	return intervals.size();
}


std::set<oneInterval*> coveredIntervals::intervals_covering_point(std::string regionID, int position) const
{
	std::set<oneInterval*> forReturn;

	if(knowRegionID(regionID))
	{
		for(std::map<int, oneInterval* >::const_iterator intervalIt = intervals_start.at(regionID).begin(); intervalIt != intervals_start.at(regionID).end(); intervalIt++)
		{
			oneInterval* existingInterval = intervalIt->second;

			if((existingInterval->regionID == regionID) && (existingInterval->from <= position) && (position <= existingInterval->to))
			{
				forReturn.insert(existingInterval);
			}
		}
	}


	return forReturn;
}

bool coveredIntervals::knowRegionID(std::string regionID) const
{
	return (bool)intervals_start.count(regionID);
}

coveredIntervals::~coveredIntervals()
{
	for(std::set<oneInterval*>::iterator iIt = intervals.begin(); iIt != intervals.end(); iIt++)
	{
		delete(*iIt);
	}

	for(std::map<std::string, IntervalTree<oneInterval*>*>::iterator treeIt = intervalTrees.begin(); treeIt != intervalTrees.end(); treeIt++)
	{
		IntervalTree<oneInterval*>* tree = treeIt->second;
		delete(tree);
	}
}

void coveredIntervals::reduceIntervalsBy(int reduction)
{
	intervalTree_upToDate = false;

	size_t intervals_before = intervals.size();
	std::set<oneInterval*> intervals_to_delete;
	for(std::set<oneInterval*>::iterator intervalIt = intervals.begin(); intervalIt != intervals.end(); intervalIt++)
	{
		oneInterval* i = *intervalIt;
		assert(i->from <= i->to);
		i->from = i->from + reduction;
		i->to = i->to - reduction;
		if(!(i->from <= i->to))
		{
			intervals_to_delete.insert(i);
		}
	}

	for(std::set<oneInterval*>::iterator intervalIt = intervals_to_delete.begin(); intervalIt != intervals_to_delete.end(); intervalIt++)
	{
		oneInterval* i = *intervalIt;
		intervals.erase(i);
		delete(i);
	}

	size_t intervals_after = intervals.size();

	intervals_start.clear();
	intervals_stop.clear();

	for(std::set<oneInterval*>::iterator intervalIt = intervals.begin(); intervalIt != intervals.end(); intervalIt++)
	{
		oneInterval* i = *intervalIt;

		assert(intervals_start[i->regionID].count(i->from) == 0);
		intervals_start[i->regionID][i->from] = i;

		assert(intervals_stop[i->regionID].count(i->to) == 0);
		intervals_stop[i->regionID][i->to] = i;
	}

	validate();

	std::cout << "coveredIntervals::reduceIntervalsBy " << reduction << ": " << intervals_before << " before, " << intervals_after << " after.\n" << std::flush;
}

void coveredIntervals::validate()
{
	std::set<oneInterval*> found_intervals_in_start_indices;
	std::set<oneInterval*> found_intervals_in_stop_indices;

	for(std::map<std::string, std::map<int, oneInterval* > >::iterator regionIt = intervals_start.begin(); regionIt != intervals_start.end(); regionIt++)
	{
		std::string regionID = regionIt->first;
		for(std::map<int, oneInterval* >::iterator positionIt = regionIt->second.begin(); positionIt != regionIt->second.end(); positionIt++)
		{
			int position = positionIt->first;
			oneInterval* i = positionIt->second;

			found_intervals_in_start_indices.insert(i);
			assert(intervals.count(i));

			assert(i->regionID == regionID);
			assert(i->from == position);
		}
	}

	for(std::map<std::string, std::map<int, oneInterval* > >::iterator regionIt = intervals_stop.begin(); regionIt != intervals_stop.end(); regionIt++)
	{
		std::string regionID = regionIt->first;
		for(std::map<int, oneInterval* >::iterator positionIt = regionIt->second.begin(); positionIt != regionIt->second.end(); positionIt++)
		{
			int position = positionIt->first;
			oneInterval* i = positionIt->second;

			found_intervals_in_stop_indices.insert(i);
			assert(intervals.count(i));

			assert(i->regionID == regionID);
			assert(i->to == position);
		}
	}

	assert(intervals.size() == found_intervals_in_start_indices.size());
	assert(intervals.size() == found_intervals_in_stop_indices.size());
	for(std::set<oneInterval*>::iterator iIt = intervals.begin(); iIt != intervals.end(); iIt++)
	{
		oneInterval* i = *iIt;
		assert(found_intervals_in_start_indices.count(i));
		assert(found_intervals_in_stop_indices.count(i));

		assert(i->from <= i->to);
	}
}
} /* namespace GraphAlignerUnique */
