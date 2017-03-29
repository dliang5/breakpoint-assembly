#include <iostream> 
#include <string>
#include <fstream>
#include <istream> 
#include <sstream> 
#include <vector> 
#include <ostream> 
#include <stdlib.h>
#include <cstdlib> 
#include <algorithm> 
#include <set> 
#include <map> 

using namespace std; 
/*
For: <primary> 65 129 <secondary> 2161 2225 
Rev: <primary> 113 177 <secondary> 2113 2177
*/

class location{ //there's sixteen columns, only use four at the end.
    public:
    string ID;
    int samFlag;
    string chrom_pos;
    int position1; //about 1000-2000 is a good position 
    int map_quality; //<trash> at the moment
    string cigarPos; //as long as it's above 20M (77) then it's all good \ 
                //if im really lazy, just convet M to int as well and subtract the size based on it.
    string equiva; //has to be '=' 
    int position2; //or the position 
    int opos; //<trash> maybe
    string nt; //the random long nucleotide sequence <trash> 
    string ntQuality; //the apparent quality checker for the sequence <trash> 
    string uk1, uk2, uk3, uk4; //unknown parts that is prob trash during the time of evaluations <trash> 
//    int counter = 0; //to see if there is atleast 3 counters then go.
    

    bool reverse ; 
    bool forward ; 

    bool operator==(const location&) const ; 
    bool operator<(const location&) const ; 

};

istream& operator>>(istream& iss, location & data){ //change this to friend operator later on and put private instead of public there
    //if(getline(iss, data.ID, '\t')){ //data ID is already given
        iss >> data.ID >> data.samFlag >> data.chrom_pos >> data.position1  
        >> data.map_quality >> data.cigarPos >> data.equiva >> data.position2 
        >> data.opos >> data.nt >> data.ntQuality >> data.uk1 >> data.uk2 
        >> data.uk3 >> data.uk4;

        return iss ; 
//        ++data.counter; //just to initialize it to 1;
    
}
//do another one for outfile operation 
//trying to find and erase common elements here
bool location::operator==(const location& rhs) const { 
    if(position1 != rhs.position1 ){
        return false ; 
    }else if(position2 != rhs.position2){
        return false ; 
    }
    return true ; 
}

bool location::operator<(const location& rhs) const { 
    return position1 < rhs.position1 ; 
}

//helper/main functions 
void algorith_search(vector< vector<location> > &, location, bool, int) ; //this will print the clusters to the files
void printCluster(vector< vector<location> > &, ofstream& ) ; //this prints the cluster to the file 
void search_potential_breaks(vector< vector<location> > &, vector< vector<location> > &, int, ofstream&, ofstream&) ; //this will use both forward and reverse clusters and try to find a potential breakpoint  
bool sortByLocation(const location &lhs , const location &rhs) ; //this sorts by position
bool compByLocation(const location &lhs , const location &rhs) ; //this compares and give me the unique location 
vector<location> duplicate( vector<location> &, vector<location> & ) ; 
void merging_clusters( vector< vector<location> > &) ; //this is just condensing the clusters into one big cluster to avoid duplications in the goodFile
vector<location> intersection(vector<location>  &, vector<location> &) ; 
vector<location> unique_difference(vector<location> & , vector<location> &) ; 
bool compByLocation1(const location &lhs , const location &rhs) ; //this compares and give me the unique location 
                //test_2 = unique_difference(merged[0], current[location_counter]) ; 
//test running program on 303_parsed.sam
bool compFirst(vector<location>, vector<location>) ; 
int main(int argc, char * argv[]){

    int cluster_size = 3000  ; 
    int min_map_quality = 20 ;
    //change the argc to 8: input output for_flag1 for_flag2 rev_flag1 rev_flag2 chrom_region 
    if(argc != 3 ){
        // printf("Usage: %s <input file> <output file>") ; 
        //modification for the current project
        cout << "Error(): Usage: ./reading <input file> <output file> " << endl ; 
        exit(1) ;
    }

    //usage example: g++ reading.cpp -o reading
    //reading 303_parsed.sam 303_sam_result.txt
    string file1 = "good_" , file2 = "garbage_" , file3 = "summary-";
    file1 += argv[2] ; file2 += argv[2] ; file3 += argv[2] ;  
    ifstream file(argv[1]); 
    ofstream outFile(argv[2], ios_base::app); //this just has all forward and reverse clusters 
    //test the current file for right now
    ofstream goodFile(file1.c_str(), ios_base::app) ; //this should have the potential breakpoints
    ofstream garbageFile(file2.c_str(), ios_base::app) ;  //everything else here and nothing important
    ofstream summary(file3.c_str(), ios_base::app) ; //this is the summary of the clusters avaliable

    int for_1 = 65 , for_2 = 129 , rev_1 = 113 , rev_2 = 177 , sec_for_1 = 2161 , sec_for_2 = 2225 , sec_rev_1 = 2113 , sec_rev_2 = 2177 ; 
    string chrom_region = "2L" ; 
    string chrom_region_1 = "X" ; 
    string chrom_region_2 = "2R" ; 
    string chrom_region_3 = "3R" ; 

    if(!file){
        cout << "The <input file> is currently empty or does not exist, please try again" << endl;
        exit(1) ; 
    }
    string line;
    vector< vector<location> > for_clusters ; //default here just to be safe for the program to work
    vector< vector<location> > rev_clusters ; 
    vector< vector<location> > potential ; 
    while(getline(file, line)){ // \n is the default delimiter here, so this goes down the input
        
        stringstream ss(line);
        location data1;
        ss >> data1;
        data1.reverse = false ; 
        data1.forward = false ; 

        //getting past headers in sam file:
        if(data1.ID.compare("@SQ") == 0){ //skipping all the header top ranging from line 1 - ~1500
            continue;
        }
 
        if ( data1.map_quality < min_map_quality || data1.equiva.compare("=") != 0) { //has to be above 20 map quality and "=" 
            continue ; //skipping the ones with bad qualities
        }
        //commented out part: uuuh reverse part will go from 3000 different reads to 10000 reads.
        if( abs(data1.position1 - data1.position2) <= 1000000 || (data1.position2 - data1.position1) < 0 ){
            continue; //skipping the ones that aren't sutiable candidates'
        }
        // cout << data1.samFlag << " this is a test" << endl ; //everything is being read fine

        //remember to add on the secondary alignments once you're done with the core purpose'
        //checking for reverse if it's there and if there are forward strands then continue leaving no other strands.
        if ( (data1.samFlag == rev_1 || data1.samFlag == rev_2 || data1.samFlag == sec_rev_1 || data1.samFlag == sec_rev_2) \
        && (data1.samFlag != for_1 || data1.samFlag != for_2 || data1.samFlag != sec_for_1 || data1.samFlag != sec_for_2) ) { 
            //this works
            data1.reverse = true ;
            data1.forward = false ;  
        }else if( (data1.samFlag == for_1 || data1.samFlag == for_2 || data1.samFlag == sec_for_1 || data1.samFlag == sec_for_2) \
        && (data1.samFlag != rev_1 || data1.samFlag !=rev_2 || data1.samFlag != sec_rev_1 || data1.samFlag != sec_rev_2) ){
            data1.reverse = false ;
            data1.forward = true ; 
        }
        else if ( (data1.samFlag != for_1 || data1.samFlag != for_2 || data1.samFlag != sec_for_1 || data1.samFlag != sec_for_2) \
         &&  (data1.samFlag != rev_1 || data1.samFlag != rev_2 || data1.samFlag != sec_rev_1 || data1.samFlag != sec_rev_2) ) { 

            garbageFile << data1.ID << "\t" << data1.samFlag << "\t" << data1.chrom_pos << "\t" << data1.map_quality << \
            "\t" << data1.position1 << "\t" << data1.position2 << "\n" ;             
            continue ; 
        }
        // this checks for the chromosomal regions that I need to check for this particular project.
        if( (data1.chrom_pos == chrom_region_2) || (data1.chrom_pos == chrom_region_1) || (data1.chrom_pos == chrom_region_2) || (data1.chrom_pos == chrom_region_3) ){

            //don't remove the comment lines until you have the merging part done with I22 or 319 whatever to make it work. '
            //do nothing here since I'm looking for these areas of the chromosomes

        }else{
            continue ; 
        }
        

        //this is the core of the function here
        //check why the forward clusters has a shit ton
        bool in_cluster = false ; 
        if( (data1.reverse == true) && (data1.forward != true) ) { 
            algorith_search(rev_clusters, data1, in_cluster,cluster_size) ;
        }else if( (data1.reverse != true) && (data1.forward == true) ){ 
            algorith_search(for_clusters, data1, in_cluster,cluster_size) ;
            //printCluster(for_clusters, garbageFile) ; 
        }

        //creating a vector that holds all of the forward first entry position to be tested with the reverse clusters later on.
        if( in_cluster == false ){
            vector<location> new_clusters ; 
            new_clusters.push_back( data1 ) ; 
            if(data1.reverse){
                rev_clusters.push_back( new_clusters ) ; 
            }else if(!data1.reverse){
                // cout << "inserting forward data here samflag " << data1.samFlag << endl ; 
                for_clusters.push_back( new_clusters ) ; 
                //potential.push_back( new_clusters ) ;
            }
            new_clusters.clear() ; 
            //cout << data1.samFlag << " this is the random spot to be here right now " << endl ; 
        } 
    }
    //merging the clusters here to make it more user readable.
    vector< vector<location> > new_for ; 
    vector< vector<location> > new_rev ;

    // printCluster(for_clusters, outFile) ; //unsensored printing here 
    // // outFile << "this is the divided line between forward and reverse clusters kldnglkandflkgndslkgkndlknglkdsnglkdnflkdsnkndslakgndslkgndlkgnadlfkngldkfnlndfdlkagnf" << endl ; 
    // outFile << "" << endl ; //for user readability here 

    // printCluster(rev_clusters, outFile) ; 
    cout << "yes" << endl ; 
    merging_clusters(for_clusters) ;
    merging_clusters(rev_clusters) ; 
    printCluster(for_clusters, outFile) ; 
    printCluster(rev_clusters, outFile) ; 
    search_potential_breaks(for_clusters, rev_clusters, cluster_size, goodFile, summary) ; 

           
//There are some files that are already separated it in a smaller scale
    file.close();
    garbageFile.close(); 
    outFile.close();
    goodFile.close(); 
    summary.close();
    return 0;
}

//functions ; 
//prints the lcuster into out file 
void printCluster(vector< vector<location> > &current, ofstream& out){ 
    for(int c = 0 ; c < current.size() ; c++){
        if(current[c].size() > 4){ 
            for(int i = 0 ; i < current[c].size() ; i++){
                out << c << "\t";
                out << current[c][i].ID << "\t" << current[c][i].samFlag << "\t" << current[c][i].chrom_pos << "\t" << current[c][i].map_quality << "\t" << current[c][i].position1 << "\t" << current[c][i].position2 << '\n';

            }
        }
    }
}

//this inputs each sequence into a cluster 
void algorith_search(vector< vector<location> > &current, location data, bool in_cluster, int clut_size){
    //means it's a reverse read 
    for (int c = 0 ; c < current.size() ; c ++ ) {
        if( current[c][0].chrom_pos != data.chrom_pos){
            continue ; 
        }else { 
            if( data.position1 < current[c][0].position1 + clut_size ){
                if ( ( abs( data.position2 - current[c][0].position2 ) ) < clut_size && ( abs( data.position1 - current[c][0].position1) ) < clut_size ){
                    current[c].push_back( data ) ; 
                    in_cluster = true ; 
                }
            }
        }
    }
}
//compare and intersect two vectors into vector c3 
//fix intersection sort here man. 
vector<location> intersection(vector<location> &c1 , vector<location> &c2){
    vector<location> c3 ; //this will be returned ; 
    sort(c1.begin(), c1.end(), sortByLocation) ;
    sort(c2.begin(), c2.end(), sortByLocation) ;
    //intersection begins here by iterating from c1 and c2 to the end by position1
    set_intersection(c1.begin(), c1.end(), 
                    c2.begin(), c2.end(), 
                    back_inserter(c3), compByLocation) ; 
    sort(c3.begin(), c3.end(), sortByLocation) ; 
    return c3 ; 
}

//compare and take out the unique elements of both vectors into vector c3
vector<location> unique_difference(vector<location> &c1, vector<location> &c2){
    vector<location> c3 ; 
    // cout <<"in here in unique difference" << endl ;
    sort(c1.begin(), c1.end(), sortByLocation) ;
    sort(c2.begin(), c2.end(), sortByLocation) ;
    //same idea as intersection but instead of looking for matching positions1, I look for non matching ones. 
    set_symmetric_difference(c1.begin(), c1.end(), 
                            c2.begin(), c2.end(), 
                            back_inserter(c3), compByLocation1) ; 
    sort(c3.begin(), c3.end(), sortByLocation) ; 
    return c3 ; 
    
}


vector<location> duplicate ( vector<location> &v1 , vector<location> &v2 ) { 
    map<location,int> read_counts ; 
    vector<location> return_vec ; 
    int count = 0 ; 
    for ( int i = 0 ; i < v1.size() ; i ++ ) { 
        read_counts[v1[i]] = 1  ; 
    }
    for ( int i = 0 ; i < v2.size() ; i ++ ) { 
        read_counts[v2[i]] ++ ; 
    }

    //// go through and find out how many values in hash are 2 
    //// if > 50% return a vector containing the keys - 
    // go down here else return nothing.
    int total = 0 ; 
    for(map<location,int>::iterator it = read_counts.begin() ; it != read_counts.end() ; ++it) {
        if(it->second == 2){ 
            count++ ;  
        }
        total ++ ; 
    }
    // determing count overlaps 50% of each vector, v1 and v2
//    if( (v1.size()) / 2 <= count && (v2.size()) / 2 <= count ) {
    if ( total/2 < count){ 
        // vector<location> return_vec ; 
        for ( map<location, int>::iterator it = read_counts.begin() ; it != read_counts.end() ; ++it ) { 
            return_vec.push_back( it->first ) ; 
        }
    }else{
        // returning a vector of 0 
        return return_vec ; 
    }
    return return_vec ; 
}


//merging the clusters here with intersection and set_difference properties of vector - is it slow? yeah.....
//it will create the first vector/cluster hold all of the same clusters as current[0] and storing the locaitons of the unmergeable clusters in an array 
void merging_clusters( vector< vector<location> > &current) {
    int erase_flag = 0 ; 
    int delete_counter = 0 ; 
    //// compare pairwise across all clusters
    for ( int i = current.size() ; i > 0 ; i -- ) { 
        for ( int i2 = i - 1 ; i2 > -1 ; i2 -- ) { 
            vector<location> dup = duplicate( current[i], current[i2] ) ;
            cout << " beforehand segfault? " << endl ;  
            if ( dup.size() > 0 ) { 
                // only go through here if it's not at the start of the program. '
                if(delete_counter != 0 && i != current.size()){
                    current.erase(current.begin() + delete_counter) ; 
                }
                current[i2] = dup ; 
                erase_flag = 1 ; 
                //segfault here  
                break ; 
            }else{ 
                erase_flag = 0 ;    
            }
        }
        if(erase_flag == 1) { 
            //saving the location to delete ; 
            delete_counter = i ; 
        }
    }
    //erasing the last one just in case 
    current.erase(current.begin() + delete_counter) ; 


/*
    vector<int> locations ; //unmergeable locations here
    vector<int> finished ;  
    vector< vector<location> > merged ; 
    //set<location> m_cur ; 
    bool has_array = false ; 
    int counter = 0 ;  
    int bookmark = 0 ; 
    while(current.size() != 0){
        //keep going until the size is 0 or something 
        // if(current[0].size() < 3){
        //     current.erase(current.begin()) ; 
        //     continue ; 
        // }
        merged.push_back(current[0]) ; 
        current.erase(current.begin()) ; //erasing the very first cluster each time 
        /// looping through the remaining vector<vector<location>> current and merge into merged[counter]
        for(int i = 0 ; i < current.size() ; i++){
            // if(current[i].size() < 3) { 
            //     finished.push_back(i) ; 
            //     continue ; 
            // }
            //checking the size of the intersection to confirm they are alike 
            vector<location> test_1, test_2 ; 
            test_1 = intersection(merged[counter], current[i]) ;
            int onefifth = merged[counter].size() / 5 ;  
            if( test_1.size() >  onefifth && test_1.size() != onefifth ){
                test_2 = unique_difference(merged[counter], current[i]) ; 
                merged[counter].insert(merged[counter].end(), current[i].begin(), current[i].end() ) ; 
                sort(merged[counter].begin(), merged[counter].end(), sortByLocation) ; 
                finished.push_back(i) ; 
            }
        }
        //removing clusters that have merged to avoid redundant clusters
        for(int j = finished.size() ; j > 0 ; j--){
            for(int i = current.size()  ; i > 0 ; i--){
                //going backwards to avoid the vector reorganizing itself based on the memory address locations 
                current.erase(current.begin() + finished[j]) ; 
        }
        finished.clear() ; 
        cout << counter << " this is the value of counter here before incrementing it" << endl ; 
        cout << "finished should be empty here " << finished.size() << endl ;
        counter++ ; 
    }
    //making sure all the duplicates are really gone here
    vector<vector<location> > now ; 
    now.resize(merged.size()) ; 
    for(int i = 0 ; i < merged.size() ; i++){
        unsigned size = merged[i].size() ; 
        set<location> m_cur(merged[i].begin(), merged[i].end()) ; 
        now[i].assign(m_cur.begin(), m_cur.end()) ; 
        m_cur.clear() ; 
    }

    // for(int i = 0 ; i < merged.size() ; i++){
    //     merged[i].erase( unique( merged[i].begin(), merged[i].end()), merged[i].end()) ; 
    // }


    return now ; 
*/

}

//I want to compare both the forward and reverse cluster and if they both are within similar distance 
//then the whole cluster is printed into the good_file.
void search_potential_breaks(vector< vector<location> > &forward , vector< vector<location> > &reverse, int clut_size, ofstream& out, ofstream& summary) {
    //use forward cluster to check reverse cluster.
    //the assumption is that the clusters are already above size 4 each.
    int counter_check = 0 ;  
    bool check = false ; 
    int move_on = 0 ; 
    for ( int f = 0 ; f < forward.size() ; f++ ){
        int flag = 0 ;
        if(forward[f].size() > 4) {
            int size = forward[f].size() ;
            int f_strand_1 = forward[f][0].position1 ;
            int f_strand_2 = forward[f][0].position2 ;
            //just to get the distance for each forward cluster and compare it to the reverse clusters and see what's up 
            for( int r = 0 ; r < reverse.size() ; r++ ){  
                if(counter_check > r){ r = counter_check ; } //skipping the ones we have as clusters . 
                int r_strand_1 = reverse[r][0].position1 ; int r_strand_2 = reverse[r][0].position2 ; 

                if (reverse[r].size() > 4) { 

                    if( ( abs( f_strand_1 - r_strand_1 ) < 30000 && abs( f_strand_1 - r_strand_1 ) < 30000  ) \
                    || ( abs( f_strand_2 - r_strand_2 ) < 30000 && abs( f_strand_2 - r_strand_2) < 30000 ) ) {
                        
                        //that might not work in general.
                        //this in theory should print out forward and reverse clusters together to show they are a good match or not
                            for(int j = 0 ; j < forward[f].size() ; j++){
                                if(j == 0 || j == forward[f].size()-1){ //for summary file
                                    summary << f << "\t";
                                    summary << forward[f][j].ID << "\t" << forward[f][j].samFlag << "\t" << forward[f][j].chrom_pos << "\t" << forward[f][j].map_quality << "\t" << forward[f][j].position1 << "\t" << forward[f][j].position2 << '\n'; 
                                }
                                out << f << "\t";
                                out << forward[f][j].ID << "\t" << forward[f][j].samFlag << "\t" << forward[f][j].chrom_pos << "\t" << forward[f][j].map_quality << "\t" << forward[f][j].position1 << "\t" << forward[f][j].position2 << '\n'; 
                            }
                            for(int j = 0 ; j < reverse[r].size() ; j++){
                                if(j == 0 || j == forward[f].size()-1){ //for summary file 
                                    summary << r << "\t";
                                    summary << reverse[r][j].ID << "\t" << reverse[r][j].samFlag << "\t" << reverse[r][j].chrom_pos << "\t" << reverse[r][j].map_quality << "\t" << reverse[r][j].position1 << "\t" << reverse[r][j].position2 << '\n'; 
                                }
                                summary << endl ; // making sure they are separated by an empty line.
                                out << r << "\t";
                                out << reverse[r][j].ID << "\t" << reverse[r][j].samFlag << "\t" << reverse[r][j].chrom_pos << "\t" << reverse[r][j].map_quality << "\t" << reverse[r][j].position1 << "\t" << reverse[r][j].position2 << '\n'; 
                            }
                            //an easier way of reading when the clusters split
                           // out << "---------------------------------------" << f << "----------------------------------------------" <<  endl ; 
                            out << endl ; 
                            counter_check = r ; 
                            break ; 
                    }
                }
            }
        }
    }
}

bool sortByLocation(const location &lhs , const location &rhs){
    if (lhs.position1 < rhs.position1){
        return true ; 
    }else if(lhs.position1 > rhs.position1){
        return false ; 
    }else if(lhs.position1 == rhs.position1){
        return false ;
    }
    // if(lhs.position2 < rhs.position2){
    //     return lhs.position2 < rhs.position2 ; 
    // }else{
    //     return lhs.position1 < rhs.position1 ;
    // } 
}
//for intersection
bool compByLocation(const location &lhs , const location &rhs){
    if(lhs.position1 == rhs.position1 && lhs.position2 == rhs.position2){
        // cout << "in true for compbylocation" << endl ; 
        return true ; 
    }else { 
        // cout << "in false for compybylocation" << endl ; 
        return false ; 
    }
}
//for unique difference
bool compByLocation1(const location &lhs , const location &rhs){
    if(lhs.position1 == rhs.position1){
        return false  ; 
    }else { 
        return true ; 
    }
}

//comparing the first positions and last positions 
bool compFirst(vector<location> c1, vector<location> c2) { 

}
