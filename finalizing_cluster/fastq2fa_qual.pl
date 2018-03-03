use strict ; 
use warnings ; 

while (<STDIN>) { 
	my $read = <STDIN> ; 
	my $trash = <STDIN> ; 
	my $qual = <STDIN> ; 

	$_ =~ s/^\@/>/ ;
	print $_, $read ; 
	
	print STDERR $_ ; 
	chomp $qual ; 
	for ( my $i = 0 ; $i < length( $qual ) ; $i ++ ) { 
		print STDERR ord( substr( $qual, $i, 1 ) ) - 33, " " ; 
	}
	print STDERR "\n" ;
}
