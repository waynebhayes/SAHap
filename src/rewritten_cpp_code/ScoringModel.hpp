#ifndef SAHAP_SCORINGMODEL_HPP
#define SAHAP_SCORINGMODEL_HPP

namespace SAHap {

/*
 * A ScoringModel provides a score for a given solution.
 * It gets invoked whenever the Genome changes to automatically keep track of the score.
 * Multiple ScoringSystems may be attached to a Genome at once.
 * 
 * Two ScoringSystems are currently implemented:
 * - WMEC:        The Weighted Minimum Error Correction model
 * - GroundTruth: Difference from the ground truth
 */
class ScoringModel {
public:
};

}

#endif