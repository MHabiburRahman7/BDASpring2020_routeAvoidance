//https://github.com/CallmeNezha/SimpleDBSCAN

#include "dbscan_nezha.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

//#define FULL_DATA
#define QUICK_TEST

using namespace std;

template <typename T>
void readCSV(string filename, vector<vector<T>>& perTrajectory, vector<T> & mixedTrajectory) {
    FILE* realfile;
    realfile = fopen(filename.c_str(), "r");

    vector<T> onelineV;

    int x;
    float lineFloat, line2Float;

#if defined(FULL_DATA)
        for (int i = 0; (x = getc(realfile)) != EOF; i++) {
            for (int j = 0; (x = getc(realfile)) != '\n'; j++) {
                fscanf(realfile, "%f,%f;", &lineFloat, &line2Float);

                //1st line float always skipped by 1 character -> and all of 1st line start with 39
                lineFloat += 39.0f;

                onelineV.push_back(T{ lineFloat, line2Float });
                mixedTrajectory.push_back(T{ lineFloat, line2Float });
            }
            perTrajectory.push_back(onelineV);
            onelineV = {};
        }
#elif defined(QUICK_TEST)
        for (int i = 0; i < 1000; i++) { //due to time remaining, so we just test this using small number of data
            for (int j = 0; (x = getc(realfile)) != '\n'; j++) {
                fscanf(realfile, "%f,%f;", &lineFloat, &line2Float);

                //1st line float always skipped by 1 character -> and all of 1st line start with 39
                lineFloat += 30.0f;

                onelineV.push_back(T{ lineFloat, line2Float });
                mixedTrajectory.push_back(T{ lineFloat, line2Float });
            }
            perTrajectory.push_back(onelineV);
            onelineV = {};
        }
#endif

    fclose(realfile);
}

template <typename T>
void writeHotspotOutput(string filename, vector<T> mixedData, vector<vector<uint>> clusterKey) {

    FILE* pFile;
    pFile = fopen(filename.c_str(), "w");

    //header for information
    fprintf(pFile, "latitude;longitude;class\n");

    for (int i = 0; i < clusterKey.size(); i++) {
        for (int j = 0; j < clusterKey[i].size(); j++) {
            fprintf(pFile, "%f;%f;%d\n", mixedData[clusterKey[i][j]][0], mixedData[clusterKey[i][j]][1], i);
        }
    }

    fclose(pFile);

}

template<typename T>
void writeTrajectoryStartEndCluster(string filename, vector<vector<uint>> the_cluster, vector<vector<T>> trajectoryData) {
    FILE* pFile;
    pFile = fopen(filename.c_str(), "w");

    //header for information
    fprintf(pFile, "trajectories; cluster_id;trajectory_id;\n");

    for (int i = 0; i < the_cluster.size(); i++) {
        for (int k = 0; k < the_cluster[i].size(); k++) {
            for (int j = 0; j < trajectoryData[the_cluster[i][k]].size(); j++) {
                fprintf(pFile, "%f;%f;", trajectoryData[the_cluster[i][k]][j][0], trajectoryData[the_cluster[i][k]][j][1]);
            }
            fprintf(pFile, "%d;%d;\n", i, the_cluster[i][k]);
        }
    }

    fclose(pFile);
}

template <typename T>
void removeOnlyOne(vector<vector<T>>& TrajectoryCollcetion, vector<int>& theScore) {

    auto it = TrajectoryCollcetion.begin();
    auto it1 = theScore.begin();

    for(; it != TrajectoryCollcetion.end() ; it++, it1++){
        if (it->size() <= 1)
        {
            TrajectoryCollcetion.erase(it--);
            theScore.erase(it1--);
        }
    }
}

template<typename T>
void catchCorrespondingData(vector<vector<T>> original_data, vector<T>& start_output, vector<T>& end_output) {
    
    for (int i = 0; i < original_data.size(); i++) {
        start_output.push_back(original_data[i][0]);
        end_output.push_back(original_data[i][original_data[i].size()-1]);
    }
}

void updateDict(vector<pair<int, int>>& the_dict, vector<vector<uint>> the_cluster, int mode) {
    //mode 0 = clustered data is start position
    //mode 1 = clustered data is end position

    for (int i = 0; i < the_cluster.size(); i++) {
        for (int j = 0; j < the_cluster[i].size(); j++) {
            if(mode == 0)
                the_dict[the_cluster[i][j]].first = i;
            else if(mode == 1)
                the_dict[the_cluster[i][j]].second = i;
        }
    }
}

template <typename T>
void writeResult(string filename,  vector<T> resultData) {
    FILE* pFile;
    pFile = fopen(filename.c_str(), "w");

    //header for information
    fprintf(pFile, "latitude;longitude;\n");

    for (int i = 0; i < resultData.size(); i++) {
        fprintf(pFile, "%f;%f;\n", resultData[i][0], resultData[i][1]);
    }

    fclose(pFile);
}

template <typename T>
void weightTheData(vector<vector<uint>> theMixedCluster, vector<int>& theScore, vector<vector<T>> fullData) {
    
    //init the score's vector
    //create range to make weighting easier
    vector<int> trajectory_range = {};
    int last_take = 0;
    for (int i = 0; i < fullData.size(); i++) {
        theScore.push_back(0);
        trajectory_range.push_back(fullData[i].size() - 1 + last_take);
        last_take += fullData[i].size()-1;
    }

    //calculate the weight
    for (int i = 0; i < theMixedCluster.size(); i++) {
        for (int j = 0; j < theMixedCluster[i].size(); j++) {
            //to check the range or this data belongs to which trajectory
            for (int k = 0; k < trajectory_range.size() - 2; k++) {
                if(trajectory_range[k] <= theMixedCluster[i][j] && theMixedCluster[i][j] <= trajectory_range[k+1])
                    theScore[k] += theMixedCluster[i].size();
            }
        }
    }
}

bool scoreSort(pair<int, int> a, pair<int, int> b) {
    return (a.second < b.second);
}

struct latLongVecf {
    float data[2];
    float operator[](int idx) const { return data[idx]; }
};

int main() {


    /**
    * @describe: Run DBSCAN clustering alogrithm
    * @param: V {std::vector<T>} : data
    * @param: dim {unsigned int} : dimension of T (a vector-like struct)
    * @param: eps {Float} : epsilon or in other words, radian
    * @param: min {unsigned int} : minimal number of points in epsilon radian, then the point is cluster core point
    * @param: disfunc {DistanceFunc} :!!!! only used in bruteforce mode.  Distance function recall. Euclidian distance is recommanded, but you can replace it by any metric measurement function
    * @usage: Object.Run() and get the cluster and noise indices from this->Clusters & this->Noise.
    * @pitfall: If you set big eps(search range) and huge density V, then kdtree will be a bottleneck of performance
    * @pitfall: You MUST ensure the data's identicality (TVector* V) during Run(), because DBSCAN just use the reference of data passed in.
    * @TODO: customize kdtree algorithm or rewrite it ,stop further searching when minimal number which indicates cluster core point condition is satisfied
    */
    // int Run(TVector* V, const uint dim, const Float eps, const uint min, const DistanceFunc& disfunc = [](const T& t1, const T& t2)->Float { return 0; });

#pragma region HotspotDetecionRegion
    //hotspot detection
    vector<vector<latLongVecf>> theData = {};
    vector<int> theDataScore = {}, mixedDataRange = {};
    vector<latLongVecf> mixedData = {};

    //readCSV("trajectory.csv", theData, mixedData, mixedDataRange);
    readCSV("trajectory.csv", theData, mixedData);

    auto dbscan = DBSCAN<latLongVecf, float>();

    dbscan.Run(&mixedData, 2, 0.01f, 2);
    auto noise = dbscan.Noise;
    auto clusters = dbscan.Clusters;

    //weight each trajectory
    weightTheData(clusters, theDataScore, theData);

    writeHotspotOutput("hotspot_class.csv", mixedData, clusters);

#pragma endregion

#pragma region ClusteringRegion
    //clean the data that only consist 1 data
    removeOnlyOne(theData, theDataScore);
    //catch the only the start and end of each trajectory
    vector<latLongVecf> startData = {}, endData = {};
    //dictionary format = start_cluster, end_cluster
    vector<pair<int, int>> dictionary_start_end;

    //build dictionary
    for (int i = 0; i < theData.size(); i++) {
        dictionary_start_end.push_back({ -1,-1 });
    }

    catchCorrespondingData(theData, startData, endData);
    //do clustering
    dbscan.Run(&startData, 2, 0.01f, 2);
    auto start_noise = dbscan.Noise;
    auto start_clusters = dbscan.Clusters;
    //update the dictionary
    updateDict(dictionary_start_end, start_clusters, 0);
    //flush the trajectory output
    writeTrajectoryStartEndCluster("start_cluster.csv", start_clusters, theData);

    dbscan.Run(&endData, 2, 0.01f, 2);
    auto end_noise = dbscan.Noise;
    auto end_clusters = dbscan.Clusters;
    //update the dictionary
    updateDict(dictionary_start_end, end_clusters, 1);
    //flush the trajectory output
    writeTrajectoryStartEndCluster("end_cluster.csv", end_clusters, theData);
    
#pragma endregion

#pragma region ReouteRecommendation
    //test case
    //start from region 3 --> 2
    int start_query = 3, end_query = 2;

    //fetch the corresponding trajectories
    vector<pair<int,int>> fetched_score_trajectoryID;
    for (int i = 0; i < dictionary_start_end.size(); i++) {
        if (dictionary_start_end[i].first == start_query && dictionary_start_end[i].second == end_query) {
            //fetched_trajectories.push_back(i);
            fetched_score_trajectoryID.push_back({i, theDataScore[i] });
        }
    }

    //sort based on score
    sort(fetched_score_trajectoryID.begin(), fetched_score_trajectoryID.end(), scoreSort);

    //give the output
    vector<latLongVecf> theResult;
    for (int i = 0; i < theData[fetched_score_trajectoryID[0].first].size(); i++) {
        theResult.push_back(theData[fetched_score_trajectoryID[0].first][i]);
    }

    writeResult("given_result.csv",theResult);

#pragma endregion

    return 0;
}