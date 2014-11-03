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

}

bool coveredIntervals::externalIntervalCovered(std::string regionID, int from, int to)
{
	std::set<oneInterval*> S = externalIntervalCoveredBy(regionID, from, to);
	assert(S.size() <= 1);
	return (S.size() > 0);
}

std::set<oneInterval*> coveredIntervals::externalIntervalCoveredBy(std::string regionID, int from, int to)
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
	assert(intervals_start.at(i->regionID).count(i->from));
	intervals_start.at(i->regionID).erase(i->from);
	assert(intervals_stop.at(i->regionID).count(i->to));
	intervals_stop.at(i->regionID).erase(i->to);
}

void coveredIntervals::addIntervalToIndices(oneInterval* i)
{
	std::string regionID = i->regionID;
	assert((intervals_start.count(regionID) == 0) || (intervals_start.at(regionID).count(i->from) == 0));
	intervals_start[i->regionID][i->from] = i;
	assert((intervals_stop.count(i->regionID) == 0) || (intervals_stop.at(i->regionID).count(i->to) == 0));
	intervals_stop[i->regionID][i->to] = i;
}    

size_t coveredIntervals::getNumIntervals()
{
	return intervals.size();
}


std::set<oneInterval*> coveredIntervals::intervals_covering_point(std::string regionID, int position)
{
	std::set<oneInterval*> forReturn;

	if(knowRegionID(regionID))
	{
		for(std::map<int, oneInterval* >::iterator intervalIt = intervals_start.at(regionID).begin(); intervalIt != intervals_start.at(regionID).end(); intervalIt++)
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

bool coveredIntervals::knowRegionID(std::string regionID)
{
	return (bool)intervals_start.count(regionID);
}

coveredIntervals::~coveredIntervals()
{
	for(std::set<oneInterval*>::iterator iIt = intervals.begin(); iIt != intervals.end(); iIt++)
	{
		delete(*iIt);
	}
}

void coveredIntervals::validate()
{
	std::set<oneInterval*> found_intervals_in_start_indices;
	std::set<oneInterval*> found_intervals_in_stop_indices;

	for(std::map<std::string, std::map<int, oneInterval* > >::iterator regionIt = intervals_start.begin(); regionIt != intervals_stop.begin(); regionIt++)
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

	for(std::map<std::string, std::map<int, oneInterval* > >::iterator regionIt = intervals_start.begin(); regionIt != intervals_stop.begin(); regionIt++)
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

		assert(i->from >= i->to);
	}
}
} /* namespace GraphAlignerUnique */
