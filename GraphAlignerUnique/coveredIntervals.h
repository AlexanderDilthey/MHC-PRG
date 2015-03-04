/*
 * coveredIntervals.h
 *
 *  Created on: 01.10.2014
 *      Author: AlexanderDilthey
 */

#ifndef COVEREDINTERVALS_H_
#define COVEREDINTERVALS_H_

#include "../intervalTree/IntervalTree.h"

#include <string>
#include <set>
#include <map>

namespace GraphAlignerUnique {

class oneInterval {
public:
	std::string regionID;
	int from;
	int to;
};

class coveredIntervals {
protected:
	std::set<oneInterval*> intervals;
	std::map<std::string, std::map<int, oneInterval* > > intervals_start;
	std::map<std::string, std::map<int, oneInterval* > > intervals_stop;
	bool intervalTree_upToDate;
	std::map<std::string, IntervalTree<oneInterval*>*> intervalTrees;
	void buildIntervalTrees();

public:
	coveredIntervals();

	void validate();

	void addPoint(std::string regionID, int position);

	std::set<oneInterval*> intervals_covering_point(std::string regionID, int position) const;

	std::set<oneInterval*> intervals_overlapping_with_givenInterval(std::string regionID, int from, int to);

	bool isThere_overlap_with_givenInterval(std::string regionID, int from, int to);


	void removeIntervalFromIndices(oneInterval* i);
	void addIntervalToIndices(oneInterval* i);

	~coveredIntervals();

	std::set<oneInterval*> externalIntervalCoveredBy(std::string regionID, int from, int to) const;
	bool externalIntervalCovered(std::string regionID, int from, int to) const;

	bool knowRegionID(std::string regionID) const;

	size_t getNumIntervals() const;

	void reduceIntervalsBy(int reduction);


};

} /* namespace GraphAlignerUnique */
#endif /* COVEREDINTERVALS_H_ */
