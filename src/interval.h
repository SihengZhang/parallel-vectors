#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>

// Helper to access infinity easily
constexpr double INF = std::numeric_limits<double>::infinity();

enum class BoundType { Open, Closed };

class Interval {
public:
    double low, high;
    BoundType lowType, highType;
    bool isEmpty;

    Interval() : low(0), high(0), lowType(BoundType::Open), highType(BoundType::Open), isEmpty(true) {}

    Interval(double l, double h, BoundType lt, BoundType ht) 
        : low(l), high(h), lowType(lt), highType(ht), isEmpty(false) {
        
        // Sanity Check: Low cannot be greater than High
        if (l > h) {
            isEmpty = true;
            return;
        }
        
        // Infinity Enforcement: Infinity must always be Open
        if (l == -INF) lowType = BoundType::Open;
        if (h == INF)  highType = BoundType::Open;

        // Check for empty specific cases like (5, 5) or [inf, inf]
        if (l == h && (lowType == BoundType::Open || highType == BoundType::Open)) {
            isEmpty = true;
        }
    }

    // Sort order: Lower starts first. If starts equal, Closed comes before Open (-inf is handled by double compare)
    bool operator<(const Interval& other) const {
        if (low != other.low) return low < other.low;
        if (lowType == BoundType::Closed && other.lowType == BoundType::Open) return true;
        return false;
    }

    // Check if intervals touch or overlap
    bool canMergeWith(const Interval& other) const {
        if (isEmpty || other.isEmpty) return false;

        // Standard Overlap
        if (high > other.low) return true;

        // Touching boundaries: [1, 2] and (2, 3] merge?
        if (high == other.low) {
            // If the touching point is finite, we need one side to be Closed.
            // If the touching point is Infinity (unlikely for "touching"), it can't be closed.
            return (highType == BoundType::Closed || other.lowType == BoundType::Closed);
        }
        return false;
    }

    Interval merge(const Interval& other) const {
        double newLow = std::min(low, other.low);
        BoundType newLowType;
        
        // Determine Low Boundary Type
        if (newLow == -INF) {
            newLowType = BoundType::Open;
        } else if (low == other.low) {
            newLowType = (lowType == BoundType::Closed || other.lowType == BoundType::Closed) ? BoundType::Closed : BoundType::Open;
        } else {
            newLowType = (low < other.low) ? lowType : other.lowType;
        }

        double newHigh = std::max(high, other.high);
        BoundType newHighType;

        // Determine High Boundary Type
        if (newHigh == INF) {
            newHighType = BoundType::Open;
        } else if (high == other.high) {
            newHighType = (highType == BoundType::Closed || other.highType == BoundType::Closed) ? BoundType::Closed : BoundType::Open;
        } else {
            newHighType = (high > other.high) ? highType : other.highType;
        }

        return Interval(newLow, newHigh, newLowType, newHighType);
    }

    Interval intersect(const Interval& other) const {
        if (isEmpty || other.isEmpty) return Interval();

        double newLow = std::max(low, other.low);
        BoundType newLowType;
        
        if (low == other.low) {
            newLowType = (lowType == BoundType::Closed && other.lowType == BoundType::Closed) ? BoundType::Closed : BoundType::Open;
        } else {
            newLowType = (low > other.low) ? lowType : other.lowType;
        }

        double newHigh = std::min(high, other.high);
        BoundType newHighType;

        if (high == other.high) {
            newHighType = (highType == BoundType::Closed && other.highType == BoundType::Closed) ? BoundType::Closed : BoundType::Open;
        } else {
            newHighType = (high < other.high) ? highType : other.highType;
        }
        
        // If intersection resulted in -INF or INF bounds, enforce Open (handled by constructor)
        return Interval(newLow, newHigh, newLowType, newHighType);
    }

    bool strictlyEndsBefore(const Interval& other) const {
        return high < other.high;
    }

    friend std::ostream& operator<<(std::ostream& os, const Interval& iv) {
        if (iv.isEmpty) return os << "{}";

        os << (iv.lowType == BoundType::Closed ? "[" : "(");
        if (iv.low == -INF) os << "-inf";
        else os << iv.low;
        
        os << ", ";
        
        if (iv.high == INF) os << "inf";
        else os << iv.high;
        os << (iv.highType == BoundType::Closed ? "]" : ")");
        return os;
    }
};

class MultiInterval {
private:
    std::vector<Interval> intervals;

    void consolidate() {
        if (intervals.empty()) return;
        std::sort(intervals.begin(), intervals.end());

        std::vector<Interval> merged;
        merged.push_back(intervals[0]);

        for (size_t i = 1; i < intervals.size(); ++i) {
            Interval& top = merged.back();
            if (top.canMergeWith(intervals[i])) {
                top = top.merge(intervals[i]);
            } else {
                merged.push_back(intervals[i]);
            }
        }
        intervals = merged;
    }

public:
    MultiInterval() {}
    MultiInterval(Interval iv) { if (!iv.isEmpty) intervals.push_back(iv); }
    MultiInterval(std::initializer_list<Interval> list) : intervals(list) { consolidate(); }

    MultiInterval Union(const MultiInterval& other) const {
        MultiInterval result;
        result.intervals = this->intervals;
        result.intervals.insert(result.intervals.end(), other.intervals.begin(), other.intervals.end());
        result.consolidate();
        return result;
    }

    MultiInterval Intersection(const MultiInterval& other) const {
        MultiInterval result;
        size_t i = 0, j = 0;
        while (i < intervals.size() && j < other.intervals.size()) {
            Interval intersection = intervals[i].intersect(other.intervals[j]);
            if (!intersection.isEmpty) result.intervals.push_back(intersection);
            if (intervals[i].strictlyEndsBefore(other.intervals[j])) i++;
            else j++;
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& os, const MultiInterval& mi) {
        if (mi.intervals.empty()) return os << "Empty Set";
        for (size_t i = 0; i < mi.intervals.size(); ++i) {
            os << mi.intervals[i];
            if (i < mi.intervals.size() - 1) os << " U ";
        }
        return os;
    }
};

#endif