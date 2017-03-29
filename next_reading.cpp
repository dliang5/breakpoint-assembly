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


//helper/main functions 
void algorith_search(vector< vector<location> > &, location, bool, int) ; //this will print the clusters to the files
void printCluster(vector< vector<location> > &, ofstream& ) ; //this prints the cluster to the file 
void search_potential_breaks(vector< vector<location> > &, vector< vector<location> > &, int, ofstream&, ofstream&) ; //this will use both forward and reverse clusters and try to find a potential breakpoint  
bool sortByLocation(const location &lhs , const location &rhs) ; //this sorts by position
bool compByLocation(const location &lhs , const location &rhs) ; //this compares and give me the unique location 
vector< vector<location> > merging_clusters( vector< vector<location> > ) ; //this is just condensing the clusters into one big cluster to avoid duplications in the goodFile
vector<location> intersection(vector<location>  , vector<location> ) ; 
vector<location> unique_difference(vector<location>  , vector<location> ) ; 
bool compByLocation1(const location &lhs , const location &rhs) ; //this compares and give me the unique location 
                //test_2 = unique_difference(merged[0], current[location_counter]) ; 
//test running program on 303_parsed.sam

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

        if( (data1.chrom_pos == chrom_region_2) /*|| (data1.chrom_pos == chrom_region_1) || (data1.chrom_pos == chrom_region_2) || (data1.chrom_pos == chrom_region_3) */){

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
    new_for = merging_clusters(for_clusters) ;
    new_rev = merging_clusters(rev_clusters) ; 
    printCluster(new_for, outFile) ; 
    printCluster(new_rev, outFile) ; 
    search_potential_breaks(new_for, new_rev, cluster_size, goodFile, summary) ; 

           
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
        if( data.position1 < current[c][0].position1 + clut_size ){
            if ( ( abs( data.position2 - current[c][0].position2 ) ) < clut_size && ( abs( data.position1 - current[c][0].position1) ) < clut_size ){
                current[c].push_back( data ) ; 
                in_cluster = true ; 
            }
        }
    }
}
//compare and intersect two vectors into vector c3 
//fix intersection sort here man. 
vector<location> intersection(vector<location> c1 , vector<location> c2){
    vector<location> c3 ; //this will be returned ; 
    cout << "in here in intersection c1 size is " << c1.size() << " c2 size is " << c2.size() << endl ; 
    if(c1.size() <= 1 || c2.size() <= 1){
        //to avoid segfaults 
        cout << "double checking two " << endl ; 
        //also theres no point in looking at size 1 since it's obviously not right.'
        return c3 ; 
    }
    cout << " double checking one " << endl ; 
    sort(c1.begin(), c1.end(), sortByLocation) ;
    sort(c2.begin(), c2.end(), sortByLocation) ;
    //intersection begins here by iterating from c1 and c2 to the end by position1
    cout << " whoa man " << endl ; 
    set_intersection(c1.begin(), c1.end(), 
                    c2.begin(), c2.end(), 
                    back_inserter(c3), compByLocation) ; 
    sort(c3.begin(), c3.end(), sortByLocation) ; 
    cout << "this is the size of c3 before transferring over " << c3.size() << endl ; 
    return c3 ; 
}
//compare and take out the unique elements of both vectors into vector c3
vector<location> unique_difference(vector<location> c1, vector<location> c2){
    vector<location> c3 ; 
    cout <<"in here in unique difference" << endl ;
    if(c1.size() <= 1 || c2.size() <= 1){
        //to avoid segfaults 
        //also theres no point in looking at size 1 since it's obviously:wq not right.'
        return c3 ; 
    }
    sort(c1.begin(), c1.end(), sortByLocation) ;
    sort(c2.begin(), c2.end(), sortByLocation) ;
    cout << " hi in the unique difference trying to trace my problem here." << endl ; 
    //same idea as intersection but instead of looking for matching positions1, I look for non matching ones. 
    set_symmetric_difference(c1.begin(), c1.end(), 
                            c2.begin(), c2.end(), 
                            back_inserter(c3), compByLocation1) ; 
    sort(c3.begin(), c3.end(), sortByLocation) ; 
    return c3 ; 
    
}

//merging the clusters here with intersection and set_difference properties of vector - is it slow? yeah.....
//it will create the first vector/cluster hold all of the same clusters as current[0] and storing the locaitons of the unmergeable clusters in an array 
vector< vector<location> > merging_clusters( vector< vector<location> > current) {
    vector<int> locations ; //unmergeable locations here
    vector<int> finished ;  
    vector< vector<location> > merged ; 
    bool has_array = false ; 
    int location_counter = 1 ;  

    while(!has_array){
        vector<location> test_1 , test_2 ; 
        // cout << "we have confirmation here" << endl ; 

        if(merged.empty()) { 
            merged.push_back(current[0]) ; //storing first cluster here 
            finished.push_back(0) ; //setting up the first index of finished here. 
        }else if( !(merged.empty()) && location_counter == current.size()){
            //come in here once the counter reaches the end to get out of loop 
            has_array = true; 
        }else{ 
            // cout << "merged[0] size is " << merged[0].size() << endl;
            // cout << "current[location_counter] size is " <<  current[location_counter].size() << endl ;       
            test_1 = intersection(merged[0], current[location_counter]) ; //checking out similar elements(?)
            //20% testing here 
            //cout << test_1.size() << endl ;
            if( test_1.size() > merged[0].size() / 5 ) {
                // test_2 = unique_difference(merged[0], current[location_counter]) ; 
                merged[0].insert(merged[0].end(), current[location_counter].begin(), current[location_counter].end() ) ; 
                sort(merged[0].begin(), merged[0].end(), sortByLocation) ; 
                merged[0].erase( unique( merged[0].begin(), merged[0].end()), merged[0].end()) ; 
                finished.push_back(location_counter) ; 
            }else{ 
                //if the 20% fails get location of current that failed
                locations.push_back(location_counter) ; //take in the spot not needed
            }
        }
        location_counter++ ; 
    }  

    //testing the next 2+ clusters and creating new ones 
    // assummption: first cluster is done, moving on to second one. 
    cout << "second time " << endl ; 
    location_counter = 0 ; //bookmark for cluster index on merged ; s
    int index_locations =  0; 
    //hmm instead of trying to fill another vector we shall delete the existing storage until it is zero. 
    while( locations.size() != finished.size()){
        //allocating the first element to the next new cluster and deleting it from 'locations' vector 
        cout << "zero?" << endl ;

        //safety measure There
        if(index_locations == locations.size() || index_locations == finished.size()){
            break ; //making sure it avoids segfaults 
        }

        int j = locations[index_locations] ; 
        cout << "one" << endl ;

        if( find(finished.begin(), finished.end(), j) != finished.end()){ 
            //moving on to the next index if it is in finished or the used cluster vector 
            continue ; 
        }else{ 
            //going on here 
            merged.push_back(current[j]) ;
            index_locations++ ;  
        }

        cout << "two" << endl ;
        //locations.erase(locations.begin() + 0) ;
        finished.push_back(locations[0]) ;  
        cout << "three" << endl; 
        location_counter++ ; //incrementing counter for merged only as a bookmark  
        for(int i = 0 ; i < locations.size() ; i++){
            //checks if the ith position is already in place,
            ////if it is, continue to the next iteration else go on. 
            cout << "double checking to find the seg fault if it is here?" ; 
            if(find(finished.begin(), finished.end(), locations[i]) != finished.end()){
                continue ; 
            }

            cout << " three.3" << endl ; 
            vector<location> test_1, test_2 ; 
            //since location has the index of the unmerged clusters in current. 
            cout << "three . 2 " << endl ; 
            cout << "location_counter is " << location_counter << " and merged size is " << merged[location_counter].size() << endl ; 
            cout << " testing to see if the variables are fine " << endl ;    
            test_1 = intersection(merged[location_counter], current[ locations[i] ]) ;
            //maybe it's from this one when test1 goes out of scope it deletes itself and that 
            //merged also has a poniter to the stuff test1 is pointing to. 
            cout << " three.1 " << endl ; 
            if(test_1.size() > merged[location_counter].size() / 5){
                //test_1.clear() ;
                //test_2 = unique_difference(merged[location_counter], current[locations[i]]) ; 
                merged[location_counter].insert( merged[location_counter].end(), test_2.begin(), test_2.end() ) ; 
                sort(merged[location_counter].begin(), merged[location_counter].end(), sortByLocation) ; 
                //deleting duplicated values.
                cout << "four" << endl; 
                //unique should be the one causing a double free corruption here. 
                //yeah it is, can't use unique erase on it man. '
                merged[location_counter].erase( unique( merged[location_counter].begin(), merged[location_counter].end() ), merged[location_counter].end() ) ;
                cout << "five" << endl ; 
                //the ones taken, delete from vector
                //this one is also invalidating the vector as well.
                finished.push_back(locations[i]) ; //putting the used clusters in here   
                //locations.erase(locations.begin() + i) ; //deleting at index positoin
                cout << "six" << endl ;
                //i-- ; //this should avoid segfaults like ex. i deletet at index 10 when the size is 11, I would iterate as 11 for the next loop which is no go. 
                cout << "seven" << endl ; 
            }
            cout << " eight " << endl ;  
        }
        cout << " nine " << endl ; 
    }


    cout << " ten " << endl ; 

    return merged ; 
}

//I want to compare both the forward and reverse cluster and if they both are within similar distance 
//then the whole cluster is printed into the good_file.
void search_potential_breaks(vector< vector<location> > &forward , vector< vector<location> > &reverse, int clut_size, ofstream& out, ofstream& summary) {
    //use forward cluster to check reverse cluster.
    //the assumption is that the clusters are already above size 4 each.
    int counter_check = 0 ; 
    vector<int> checkoff_max ;
    vector<int> checkoff_min ; //a more strict way to avoid duplications 
    bool check = false ; 
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
    if(lhs.position2 < rhs.position2){
        return lhs.position2 < rhs.position2 ; 
    }else{
        return lhs.position1 < rhs.position1 ;
    } 
}
//for intersection
bool compByLocation(const location &lhs , const location &rhs){
    if(lhs.position1 == rhs.position1){
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