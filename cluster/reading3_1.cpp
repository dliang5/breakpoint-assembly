/* This clusters drosophila dna sequences with the chrom, samFlag, and both poositions with a window size of 2000 bases.
   It then fixes the issue of having dupliated clusters by merging them into one for the most part. 

   compile : g++ reading3_1.cpp -o reading (no flags)
   run :  ./reading sam_name sam->Ral_name (ex. ./reading 303.sam 303-results)
   output : good_<Ral_name>-result, good_<Ral_name>-summary, <Ral_name>-results


    UNDER CONSTRUCTION, THIS WORKS FOR ONE PROGRAM NOT MULTIPLE
*/

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
    string equiva ; 
    string cigarPos; //as long as it's above 20M (77) then it's all good \ 
    string equiva; //has to be '=' 
    int position2; //or the position 
    int opos; //<trash> maybe
    string nt; //the random long nucleotide sequence <trash> 
    string ntQuality; //the apparent quality checker for the sequence <trash> 
    string uk1, uk2, uk3, uk4; //unknown parts that is prob trash during the time of evaluations <trash> 
    bool reverse ; 
    bool forward ; 

    bool operator==(const location&) const ; 
    bool operator<(const location&) const ; 

};

// this reads the line and converts everything for the class
istream& operator>>(istream& iss, location & data){ //change this to friend operator later on and put private instead of public there
    //if(getline(iss, data.ID, '\t')){ //data ID is already given
        iss >> data.ID >> data.samFlag >> data.chrom_pos >> data.position1  
        >> data.map_quality >> data.cigarPos >> data.equiva >> data.position2 
        >> data.opos >> data.nt >> data.ntQuality >> data.uk1 >> data.uk2 
        >> data.uk3 >> data.uk4;

        return iss ;     
}

bool location::operator<(const location& rhs) const { 
    return position1 < rhs.position1 ; 
}

// comparing by position1 for the smallest
bool sortByLocation(const location &lhs , const location &rhs){
    if (lhs.position1 < rhs.position1){
        return true ; 
    }else if(lhs.position1 > rhs.position1){
        return false ; 
    }else if(lhs.position1 == rhs.position1){
        return false ;
    }
}
// comparing by postion2 for the smallest
bool sortByLocation2(const location &lhs, const location &rhs){
    if (lhs.position2 < rhs.position2){
        return true ; 
    }else if (lhs.position2 > rhs.position2){
        return false ; 
    }else if(lhs.position2 == rhs.position2){
        return false ; 
    }
}

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

// reads each sequence and input them into a cluster
//// TODO: already selects what kind of chromosome the vector is down there so no need for if statement double check tho 
// if cur_sequence matches a cluster that is similar, set in_cluster = true
void cluster_search(vector< vector<location> > &current, location cur_sequence, bool in_cluster, int clut_size){
    //means it's a reverse read 
    for (int c = 0 ; c < current.size() ; c ++ ) {
        if( current[c][0].chrom_pos != cur_sequence.chrom_pos){
            continue ; 
        }else { 
            if( cur_sequence.position1 < current[c][0].position1 + clut_size ){
                if ( ( abs( cur_sequence.position2 - current[c][0].position2 ) ) < clut_size && ( abs( cur_sequence.position1 - current[c][0].position1) ) < clut_size ){
                    current[c].push_back( cur_sequence ) ; 
                    in_cluster = true ; 
                }
            }
        }
    }
}

// checks two vectors and see if they match for the most part
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

// merging similar clusters
void merging_clusters( vector< vector<location> > &current) {
    int erase_flag = 0 ; 
    int delete_counter = 0 ; 
    //// compare pairwise across all clusters
    for ( int i = current.size() ; i > 0 ; i -- ) { 
        for ( int i2 = i - 1 ; i2 > -1 ; i2 -- ) { 
            vector<location> dup = duplicate( current[i], current[i2] ) ;
            //cout << " beforehand segfault? " << endl ;  
            if ( dup.size() > 0 ) { 
                // only go through here if it's not at the start of the program. '
                if(delete_counter != 0 && i != current.size()){
                    current.erase(current.begin() + delete_counter) ; 
                }
                current[i2] = dup ; 
                erase_flag = 1 ; 
                break ; 
            }else{ 
                erase_flag = 0 ;    
            }
        }
        if(erase_flag == 1) { 
            //saving the location to delete to avoid segfaults; 
            delete_counter = i ; 
        }
    }
    //erasing the last one since the iteration doesn't go back and delete it
    current.erase(current.begin() + delete_counter) ; 
}


void ordered_cluster_print(vector<vector<location> > &forward , vector<vector<location> > &reverse, int clut_size, ofstream& file_out) {
    int counter_check = 0 ;   
    for ( int f = 0 ; f < forward.size() ; f++ ){
        int flag = 0 ;
        if(forward[f].size() > 4) {
            int size = forward[f].size() ;
            int f_strand_1 = forward[f][0].position1 ;
            int f_strand_2 = forward[f][0].position2 ;
            //just to get the distance for each forward cluster and compare it to the reverse clusters and see what's up 
            for( int r = 0 ; r < reverse.size() ; r++ ){  
                if(counter_check > r){ 
                    r = counter_check ; 
                } //skipping the ones we have as clusters . 
                
                int r_strand_1 = reverse[r][0].position1 ; int r_strand_2 = reverse[r][0].position2 ; 

                if (reverse[r].size() > 4) {
                    if( ( abs( f_strand_1 - r_strand_1 ) < 30000 && abs( f_strand_1 - r_strand_1 ) < 30000  ) \
                    && ( abs( f_strand_2 - r_strand_2 ) < 30000 && abs( f_strand_2 - r_strand_2) < 30000 ) ) {

                        for(int j = 0 ; j < forward[f].size() ; j++){
                            file_out << f << "\t";
                            file_out << forward[f][j].ID << "\t" << forward[f][j].samFlag << "\t" << forward[f][j].chrom_pos << "\t" << forward[f][j].map_quality << "\t" << forward[f][j].position1 << "\t" << forward[f][j].position2 << '\n'; 
                        }
                        for(int j = 0 ; j < reverse[r].size() ; j++){
                            file_out << r << "\t";
                            file_out << reverse[r][j].ID << "\t" << reverse[r][j].samFlag << "\t" << reverse[r][j].chrom_pos << "\t" << reverse[r][j].map_quality << "\t" << reverse[r][j].position1 << "\t" << reverse[r][j].position2 << '\n'; 
                        }
                        file_out << endl ; 
                        counter_check = r ; 
                        break ; 
                    }
                }
            }
        }
    }
}

// same as up but make a summary of them
// needs its own function call here since summary has to be sorted for readability
void produce_summary(vector<vector<location> > &forward, vector<vector<location> > &reverse, ofstream& summary) { 
    int counter_check = 0 ;  
    
    for ( int f = 0 ; f < forward.size() ; f++ ){
        int flag = 0 ;

        if(forward[f].size() > 4) {
            int size = forward[f].size() ;
            int f_strand_1 = forward[f][0].position1 ;
            int f_strand_2 = forward[f][0].position2 ;
            for( int r = 0 ; r < reverse.size() ; r++ ){  

                // to avoid doing the same thing twice
                if(counter_check > r){ 
                    r = counter_check ; 
                    
                }  
                int r_strand_1 = reverse[r][0].position1 ; int r_strand_2 = reverse[r][0].position2 ; 

                // the idea is to sort both forward and reverse clusters from smallest to biggest
                if (reverse[r].size() > 4) { 

                    if( ( abs( f_strand_1 - r_strand_1 ) < 30000 && abs( f_strand_1 - r_strand_1 ) < 30000  ) \
                    && ( abs( f_strand_2 - r_strand_2 ) < 30000 && abs( f_strand_2 - r_strand_2) < 30000 ) ) {

                        
                        int f_size = forward[f].size() ;
                        int r_size = reverse[r].size() ;
                        summary << f << '\t' << forward[f][0].samFlag << '\t' << forward[f][0].chrom_pos << '\t' << forward[f][0].position1 ; 
                        int position1 = 0, position2 = 0 ; 
                        position1 = forward[f][f_size-1].position1 ; //the highest position 
                        sort(forward[f].begin(), forward[f].end(), sortByLocation2) ;
                        summary << '\t' << forward[f][0].position2 ; 
                        position2 = forward[f][f_size-1].position2 ;
                        summary << '\t' << position1 << '\t' << position2 << '\t'<< f_size;  
                        summary << endl ; 

                        summary << r << '\t' << reverse[r][0].samFlag << '\t'  << reverse[r][0].chrom_pos << '\t' << reverse[r][0].position1 ; 
                        sort(reverse[r].begin(), reverse[r].end(), sortByLocation2) ;
                        summary << '\t' << reverse[r][0].position2 ; 
                        sort(reverse[r].begin(), reverse[r].end(), sortByLocation) ; 
                        summary << '\t' << reverse[r][r_size-1].position1 ; 
                        sort(reverse[r].begin(), reverse[r].end(), sortByLocation2) ; 
                        summary << '\t' << reverse[r][r_size-1].position2 << '\t' << r_size;  

                        counter_check = r ; 
                        summary << endl << endl ; 
                        break ; 
                    }
                }
            }
        }
    }
}

int main(int argc, char * argv[]){

    int cluster_size = 3000  ; // window size 
    int min_map_quality = 20 ; // min quality for each sequence 

    if(argc != 2 ){
        cout << "Error(): Usage: ./reading <input file> " << endl ; 
        exit(1) ;
    }

    string file1 = "good_" ; 
    string file3 = "summary-";
    string args = argv[1] ; 
    if ( string::npos == args.find_first_of(".") ) {
        int period = args.find(".") ; 
        string new_file = args.substr(0, period) ; 
        file1 += new_file ; 
        file3 += new_file ; 
    }else { 
        file1 += argv[1] ; 
        file3 += argv[1] ; 
    }

    ifstream file(argv[1]); 
    //test the current file for right now
    ofstream goodFile(file1.c_str(), ios_base::app) ; //this should have the potential breakpoints
    ofstream summary(file3.c_str(), ios_base::app) ; //this is the summary of the clusters avaliable

    // samflags to check out
    int for_1 = 65 , for_2 = 129 , rev_1 = 113 , rev_2 = 177 , 
    sec_for_1 = 2161 , sec_for_2 = 2225 , sec_rev_1 = 2113 , sec_rev_2 = 2177 ; 

    string chrom [5] = {"2L", "3L", "X", "2R", "3R"} ; 

    // setting up 
    map<string, vector<vector<location> > > for_clusters ; 
    map<string, vector<vector<location> > > rev_clusters ; 

    if(!file){
        cout << "The <input file> is currently empty or does not exist, please try again" << endl;
        exit(1) ; 
    }
    string line;
    // vector< vector<location> > for_clusters ; //default here just to be safe for the program to work
    // vector< vector<location> > rev_clusters ;
    
    while(getline(file, line)){ // \n is the default delimiter here, so this goes down the input
        int counter = 0 ; 
        stringstream ss(line);
        location cur_sequence;
        ss >> cur_sequence;
        cur_sequence.reverse = false ; 
        cur_sequence.forward = false ; 

        //getting past headers in sam file:
        if(cur_sequence.ID.compare("@SQ") == 0){ //skipping all the header top ranging from line 1 - ~1500
            continue;
        }
        
        //skipping bad quality or no matching pairs
        if ( cur_sequence.map_quality < min_map_quality || cur_sequence.equiva.compare("=") != 0) {  
            continue ; 
        }

        // positions has to not be negative or less than a 1 Mb
        if( abs(cur_sequence.position1 - cur_sequence.position2) <= 1000000 || (cur_sequence.position2 - cur_sequence.position1) < 0 ){
            continue; //skipping the ones that aren't sutiable candidates'
        }

        // checking what chrom the selected sequence is
        if(cur_sequence.chrom_pos == chrom[0]){ counter = 0;} 
        else if (cur_sequence.chrom_pos == chrom[1]){ counter = 1;} 
        else if (cur_sequence.chrom_pos == chrom[2]){ counter = 2;}  
        else if (cur_sequence.chrom_pos == chrom[3]){ counter = 3;}
        else if (cur_sequence.chrom_pos == chrom[4]){ counter = 4;
        }else{           
            continue ; 
        }

        if ( (cur_sequence.samFlag == rev_1 || cur_sequence.samFlag == rev_2 || cur_sequence.samFlag == sec_rev_1 || cur_sequence.samFlag == sec_rev_2) ) { 
            cur_sequence.reverse = true ;
            cur_sequence.forward = false ;  
        }else if( (cur_sequence.samFlag == for_1 || cur_sequence.samFlag == for_2 || cur_sequence.samFlag == sec_for_1 || cur_sequence.samFlag == sec_for_2) ){
            cur_sequence.reverse = false ;
            cur_sequence.forward = true ;
            
        //// triple checking that nothing goes through here
        }else if ( (cur_sequence.samFlag != for_1 || cur_sequence.samFlag != for_2 || cur_sequence.samFlag != sec_for_1 || cur_sequence.samFlag != sec_for_2) \
         &&  (cur_sequence.samFlag != rev_1 || cur_sequence.samFlag != rev_2 || cur_sequence.samFlag != sec_rev_1 || cur_sequence.samFlag != sec_rev_2) ) {        
            continue ; 
        }

        //clustering section here
        bool in_cluster = false ; 
        if( (cur_sequence.reverse == true) && (cur_sequence.forward != true) ) { 
            cluster_search(rev_clusters[chrom[counter]], cur_sequence, in_cluster,cluster_size) ;

        }else if( (cur_sequence.reverse != true) && (cur_sequence.forward == true) ){ 
            cluster_search(for_clusters[chrom[counter]], cur_sequence, in_cluster,cluster_size) ;
        }
        
        // if a really different sequence positions come, set its own cluster.
        if( in_cluster == false ){
            vector<location> new_clusters ; 
            new_clusters.push_back( cur_sequence ) ; 
            if(cur_sequence.reverse){
                rev_clusters[chrom[counter]].push_back( new_clusters ) ; 
            }else if(!cur_sequence.reverse){
                for_clusters[chrom[counter]].push_back( new_clusters ) ;
            }
            new_clusters.clear() ; 
        }  
    }
    //merging and writing the clusters to the files.
    //trying to make this more dynamic to changes
    for (int i = 0 ; i < (sizeof(chrom) / sizeof(*chrom)) ; i++){
        merging_clusters(for_clusters[chrom[i]]) ;
        merging_clusters(rev_clusters[chrom[i]]) ; 
        ordered_cluster_print(for_clusters[chrom[i]], rev_clusters[chrom[i]], cluster_size, goodFile) ; 
        produce_summary(for_clusters[chrom[i]], rev_clusters[chrom[i]], summary) ; 
    }           
//There are some files that are already separated it in a smaller scale
    file.close();
    goodFile.close(); 
    summary.close();
    return 0;
}