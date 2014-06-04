#pragma once
#include <vector>
#include <iostream>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;
using std::max;

template <class T>
class StatisticalValue
{
private:
    vector<T> m_currentValue;
    vector<T> m_sum;
    vector<T> m_sumSquared;
    int m_numberOfSamples;
    int m_numberOfBins;
public:
    StatisticalValue(int numberOfElements = 1) :
        m_numberOfSamples(0),
        m_numberOfBins(numberOfElements)
    {
        m_currentValue.resize(numberOfElements);
        m_sum.resize(numberOfElements);
        m_sumSquared.resize(numberOfElements);
    }

    vector<T> currentValue()
    {
        return m_currentValue;
    }

    T currentValueScalar()
    {
        return m_currentValue.at(0);
    }

    vector<T> sum()
    {
        return m_sum;
    }

    vector<T> sumSquared()
    {
        return m_sumSquared;
    }

    int numberOfSamples()
    {
        return m_numberOfSamples;
    }

    int numberOfBins()
    {
        return m_numberOfBins;
    }

    void resize(int numberOfBins) {
        m_numberOfBins = numberOfBins;
        m_currentValue.resize(numberOfBins,0);
        m_sum.resize(numberOfBins,0);;
        m_sumSquared.resize(numberOfBins,0);;
    }

    vector<T> getStandardDeviation() {
        vector<T> variance = getVariance();
        vector<T> standardDeviation = variance;
        for(T &value : standardDeviation) {
            value = sqrt(value);
        }

        return standardDeviation;
    }

    vector<T> getVariance() {
        vector<T> average = getAverage();
        vector<T> squaredAverage = getSquaredAverage();
        vector<T> variance(average.size());

        for(int i=0; i<average.size(); i++) {
            variance[i] = squaredAverage[i] - average[i]*average[i];
        }

        return variance;
    }

    vector<T> getSquaredAverage() {
        vector<T> squaredAverage = m_sumSquared;
        for(T &value : squaredAverage) {
            value /= max(m_numberOfSamples, 1); // Divide by at least 1 to make sure we don't get any NaN values
        }

        return squaredAverage;
    }

    vector<T> getAverage() {
        vector<T> average = m_sum;
        for(T &value : average) {
            value /= max(m_numberOfSamples, 1); // Divide by at least 1 to make sure we don't get any NaN values
        }

        return average;
    }

    void addValue(vector<T> value) {
        m_currentValue = value;

        // Assume std vector
        for(int i=0; i<value.size(); i++) {
            m_sum[i] += value[i];
            m_sumSquared[i] += value[i]*value[i];
        }

        m_numberOfSamples++;
    }

    void addValue(T value) {
        if(m_sum.size() > 1) cout << "Warning, appending scalar to vector statistical value." << endl;
        m_currentValue.at(0) = value;
        m_sum.at(0) += value;
        m_sumSquared.at(0) += value*value;
        m_numberOfSamples++;
    }
};
