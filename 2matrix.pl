#!/usr/bin/perl 



############################### This program is free software; you can redistribute it and/or modify
############################### it under the terms of the GNU General Public License as published by
############################### the Free Software Foundation; either version 2 of the License, or
############################### (at your option) any later version.
############################### 
############################### This program is distributed in the hope that it will be useful,
############################### but WITHOUT ANY WARRANTY; without even the implied warranty of
############################### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
############################### GNU General Public License for more details.
############################### 
############################### You should have received a copy of the GNU General Public License
############################### along with this program; if not, write to the Free Software
############################### Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
###############################
###############################
###############################
############################### Copyright 2013 Nelson R. Salinas and Damon P. Little



############################### MODULES
use strict;



############################### CHECK INPUT
my $aa = ();
my $outputName = ();
my $phylip = 0;
my $nexus = 0;
my $xread = 0;
my $codeIndel = 1;
my $morpho = ();
my @infile = ();
my $totalChars = ();
my $part = ();
my $partSize = ();
my $partIndels = ();
my @stem = ();
my %stem = ();
my $debug = ();
my @taxa = ();
my @bases = ();
my $add = ();
my $RAXMLpoly = 0;

for(my $k = $#ARGV; $k >= 0; $k--){
	if($ARGV[$k] eq '-i'){
		if(-e $ARGV[$k+1]){
			push(@infile, $ARGV[$k+1]);
			} else{
				die("Could not read $ARGV[$k+1].\n");
				}
		next;		
		}
	if($ARGV[$k] eq '-n'){
		$outputName = $ARGV[$k+1];
		next;
		}
	if($ARGV[$k] eq '-d'){
		$codeIndel = 0;
		next;
		}
	if($ARGV[$k] eq '-o'){
		if($ARGV[$k+1] =~ m/x/i){
			$xread = 1;
			}
		if($ARGV[$k+1] =~ m/n/i){
			$nexus = 1;
			}
		if($ARGV[$k+1] =~ m/p/i){
			$phylip = 1;
			}
		if(!($xread || $nexus || $phylip)){
			die("You must set an output format.\n");
			}
		next;	
		}
	if($ARGV[$k] eq '-s'){
		if(length($ARGV[$k+1])){
			push(@stem, $ARGV[$k+1]);
			}
		}
	}
my $j = $#stem;
for(my $k = $#infile; $k >= 0; $k--){
	if($infile[$k] =~ m/\.fasta$|\.faa$|\.fst$|\.fa$|\.fas$|\.fna$|\.ffn$|\.frn$/){
		if($#stem == -1){
			$stem{$infile[$k]} = 'sequence';
			} elsif(length($stem[$j])){
				$stem{$infile[$k]} = $stem[$j];
				$j--;
				} else {
					die("The number of character stem names does not match the number of input files. Character names will probably not make much sense.\n");
					}
		if($infile[$k] =~ m/\.faa$/){
			$aa->[$k] = 1;
			} else {
				$aa->[$k] = 0;
				}
		} else {
			$aa->[$k] = 0;
			}
	}
undef(@stem);	


	
if(($#infile != -1) && ($xread || $nexus || $phylip) && length($outputName)){ ############################### START
	my $tTpTd = {}; ### hash of terminal => partition name => data as a string
	my $pTc = {}; ### hash of partion name => array of character names
	
	
	my($day, $month, $year) = (localtime)[3,4,5];
	my @months = ('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');
	my $quote = 'Data processed on ' . ($year + 1900) . ' ' . $months[$month] . ' ' .  $day . ' using 2matrix.pl ' . join(' ', @ARGV);
	my $quoteIndel = ();

	for(my $f = $#infile; $f >= 0; $f--){ ############################### READ INFILE(S)	
		my $matrix->[0][0] = ();
		my $terms = 0;
		open(INFILE, '<', $infile[$f]) || die("Could not open $infile[$f]!\n");
		my $firstline = <INFILE>;
		close(INFILE);
		open(INFILE, '<', $infile[$f]) || die("Could not open $infile[$f]!\n");
		if($firstline =~ m/xread/){ ############################### INPUT IN NONA FORMAT
			my $gotQuote = 0;
			my $gotChars = 0;
			my $gotMatrix = 0;
			my $quoteStart = 0;
			my $characters = ();
			my $charNames = 0;
			my @names = ();
			while(my $line = <INFILE>){
				$line =~ s/[\r\n]+$//;				
				if($gotChars == 0){
					my @letters = split('', $line);
					for(my $v = 0; $v <= $#letters; $v++){
						if($gotQuote == 0){
							if($letters[$v] eq "'"){
								if($quoteStart == 1){ ### last quote
									$gotQuote = 1;
									} else {
										$quoteStart = 1;
										}
								} elsif($quoteStart == 1){
									$quote .= $letters[$v];
									}
							} else {
								if(($letters[$v] =~ /\d/) && ($gotChars == 0)){ 
									$characters .= $letters[$v];
									}
								if(($letters[$v] eq ' ') && length($characters)){
									$gotChars = 1;
									}
								}
						}
					} elsif($gotMatrix == 0) { ### matrix lines					
						if($line =~ m/^;/){
							$gotMatrix = 1;
							} else {
								my @words = split(/ /, $line);  ### $words[0] is OTU's name
								$words[0] =~ tr/ \/\-/\_\_\_/;
								$words[0] =~ tr/[A-z][0-9]\_\.//cd;
								$matrix->[$terms][0] = $words[0];
								my $i = 1;
								my $polymorphism = 0;
								for(my $v = 1; $v <= $#words; $v++){
									my @atom = split(//, $words[$v]);
									for(my $w = 0; $w <= $#atom; $w++){
										if($atom[$w] eq '['){
											$polymorphism = 1;
											} elsif(($polymorphism == 1) && ($atom[$w] eq ']')){
												$polymorphism = 0;
												}
										$matrix->[$terms][$i] .= $atom[$w];
										if($polymorphism == 0){
											$i++;
											}
										}
									}
								$terms++;	
								}
						} elsif($line !~ m/xread/){
							if(($line =~ m/(^cc)(\s)*(\+|\-)((\s)*(\+|\-)*(\s)*(\.|\d)(\s)*)+(;$)/) && ($line !~ m/\[|\]/)){
								$line =~ tr/[0-9]\. \+\-//cd;
								$line =~ s/\+/ \+ /g;
								$line =~ s/\-/ \- /g;
								$line =~ s/\./ \. /g;
								$line =~ tr/ / /s;
								my @bits = split(/ /, $line);
								my $state = ();
								for(my $g = 0; $g <= $#bits; $g++){
									if($bits[$g] eq '-'){
										$state = 0;
										next;
										}
									if($bits[$g] eq '+'){
										$state = 1;
										next;
										}
									if($bits[$g] =~ m/\d/){
										if($bits[$g+1] eq '.'){
											my $end = $characters;
											if($bits[$g+2] =~ m/\d/){
												$end = $bits[$g+2];
												}
											for(my $h = $end; $h >= $bits[$g]; $h--){
												$add->[$f][$h] = $state;
												}
											} else {
												$add->[$f][$bits[$g]] = $state;												
												}
										next;
										}
									if(($bits[$g] eq '.') && ($bits[$g-1] !~ m/\d/)){
										my $end = $characters;
										if($bits[$g+1] =~ m/\d/){
											$end = $bits[$g+1];
											}
										for(my $h = $end; $h >= 0; $h--){
											$add->[$f][$h] = $state;
											}																			
										}
									}
								} elsif($line =~ m/^cn/){ ### cn {0 seq1 A C G T /;
									$charNames = 1;
									$line =~ s/cn {\d+ //;
									push(@names, $line);
									} elsif($charNames == 1){
										if($line =~ m/^{\d+/){
											$line =~ s/{\d+ //;
											push(@names, $line);
											} else { ### end of character/state names
												last;
												}
										}
							@{$pTc->{$infile[$f]}} = @names;
							$partSize->[$f] = $characters;
							}
				}
			for(my $h = $characters; $h >= 0; $h--){
				if($add->[$f][$h] != 0){
					$add->[$f][$h] = 1;
					}
				}
			for(my $k = $#{$matrix}; $k >= 0; $k--){
				if(exists($tTpTd->{$matrix->[$k][0]}->{$infile[$f]})){
					die("Terminals must have unique names ($matrix->[$k][0] is duplicated in $infile[$f])!\n");			
					}
				for(my $j = 1; $j <= $#{$matrix->[$k]}; $j++){
					$tTpTd->{$matrix->[$k][0]}->{$infile[$f]} .= $matrix->[$k][$j];
					}
				}
			$morpho->[$f] = 1;		
			next;
			} elsif(($firstline =~ m/,/) && ($infile[$f] =~ m/.csv$/) && ($firstline !~ m/^>/)){
				my $statesOrder = {};
				my $missingError = 0;
				chomp($firstline);
				$firstline =~ tr/ /_/;
				
				
				
				
				my @charNames = ();
				sub antiQuote {
					my $quote = $_[0];
					my $line = $_[1];
					my @output = ();
					my @bits = split(/,/, $line);
					my $quoteStart = -1;
					my $quoteEnd = -1;
					for(my $c = 0; $c <= $#bits; $c++){
						if($bits[$c] =~ m/^$quote/){
							$quoteStart = $c;							
							}
						if(($quoteStart >= 0) && ($bits[$c] =~ m/$quote$/)){
							$quoteEnd = $c;
							} 
						if((($quoteStart = -1) && ($quoteEnd = -1)) || (($quoteStart == $quoteEnd) && ($quoteStart >= 0))){
							push(@output, $bits[$c]);
							$quoteStart = -1;
							$quoteEnd = -1;
							next;
							}
						if(($quoteStart >= 0) && ($quoteEnd >= 0) && ($quoteEnd > $quoteStart)){
							push(@output, join('', @bits[$quoteStart..$quoteEnd]));	
							$quoteStart = -1;
							$quoteEnd = -1;
							} 
						}
					my $output = join(',', @output);
					$output =~ s/$quote//g;
					return($output);				
					}
				if($firstline =~ m/"/){
					@charNames = split(/,/, antiQuote('"', $firstline));						
					} else {
						@charNames = split(/,/, $firstline);						
						}
				$partSize->[$f] = $#charNames;		
				my $fieldCounter = $#charNames;		
				shift(@charNames);
				my $lineCounter = 1;
				while(my $line = <INFILE>){
					$line =~ s/[\r\n]+$//;
					$line =~ tr/\/\-/\_\_/;
					$line =~ tr/[A-z][0-9]\_\ ",?\.//cd;
					if(length($line)){
						$debug .= "Line $lineCounter: $line\n";
						my @fields = ();
						if($line =~ m/"/){
							@fields = split(/,/, antiQuote('"', $line));
							} else {
								@fields = split(/,/, $line);
								}
						if($lineCounter == 2){
							for(my $w = 1; $w <= $#fields; $w++){
								if($fields[$w] =~ m/non_additive|unordered/i){
									$add->[$f][$w-1] = 0;
									} elsif((length($fields[$w]) >= 2) && ($fields[$w] =~ m/ /)){										
										my @additiveStates = split(/ /,$fields[$w]);
										for(my $x = $#additiveStates; $x >= 0; $x--){
											$statesOrder->{$additiveStates[$x]} = $x;
											}
										$add->[$f][$w-1] = 1;
										} else {
											die("Accepted labels for character additivity are \"additive\", \"non-additive\", \"ordered\", and \"unordered\"\n\n"); 
											}
								}
							} elsif($lineCounter > 2){
								if(($fieldCounter != $#fields) && ($#fields >= 1)){
									die("It appears that rows have different number of columns in $infile[$f] (line $lineCounter).\nThis may be caused by the use of commas or double quotes within cells. Attention\nto detail is very important in systematics. If you repeatably get this error,\nyou may wish to change careers.\n\n");
									}
								$matrix->[$lineCounter-3][0] = $fields[0];	
								for(my $c = 1; $c <= $#fields; $c++){								
									$matrix->[$lineCounter-3][$c] = ' ' . $fields[$c] . ' ';
									}
								}
						$lineCounter++;
						}
					}
				for(my $d = 1; $d <= $#{$matrix->[1]}; $d++){
					my $charStates = {};
					my $stateCounter = 0;
					for(my $e = $#{$matrix}; $e >= 0; $e--){
						if($matrix->[$e][$d] =~ m/"/){
							die("It appears that rows were read with different number of columns in $infile[$f] (line $e).\nThis may be caused by the use of commas or double quotes within cells. Attention\nto detail is very important in systematics. If you repeatably get this error,\nyou may wish to change careers.\n\n");
							}
						if(($matrix->[$e][$d] eq ' ? ') || ($matrix->[$e][$d] eq ' _ ')){
							$matrix->[$e][$d] =~ tr/\_/\-/;
							$matrix->[$e][$d] =~ tr/ //d;
							next;							
							}
						if($add->[$f][$d-1] == 1){ ### character is additive
							my @tokens = split(/ /, $matrix->[$e][$d]);
							$matrix->[$e][$d] = ();
							for(my $t = $#tokens; $t >= 0; $t--){
								if(length($tokens[$t])){
									if(exists($statesOrder->{$tokens[$t]})){
										$matrix->[$e][$d] .= $statesOrder->{$tokens[$t]};
										$charStates->{$tokens[$t]} = $statesOrder->{$tokens[$t]};
										} else {
											die("Every state is sacred. State '$tokens[$t]' is not included in the additivity row for character $d.\n\n");										
											}
									}	
								}
							} else{ ### character is non-additive
								my @tokens = split(/ /, $matrix->[$e][$d]);
								for(my $t = $#tokens; $t >= 0; $t--){
									if(length($tokens[$t])){
										if(!exists($charStates->{$tokens[$t]})){
											$charStates->{$tokens[$t]} = $stateCounter;
											$stateCounter++;
											}
										$matrix->[$e][$d] =~ s/ $tokens[$t] / $charStates->{$tokens[$t]} /g;
										}
									}
								}
						$matrix->[$e][$d] =~ tr/ //d;
						if(length($matrix->[$e][$d]) >= 2){
							$matrix->[$e][$d] = '[' . $matrix->[$e][$d] . ']';
							}
						if(length($matrix->[$e][$d]) == -1){
							$matrix->[$e][$d] = '?';
							if($missingError == 0){
								print(STDERR "At least one empty cell was found. Empty cells were converted to missing ('?').\n");
								$missingError = 1;
								}
							}
						}
					my @charKeys = sort({$charStates->{$b} <=> $charStates->{$a}} keys(%{$charStates}));
					for(my $p = $#charKeys; $p >= 0; $p--){
						$charNames[$d-1] .= ' ' . $charKeys[$p];
						}
					push(@{$pTc->{$infile[$f]}}, "$charNames[$d-1];");
					}
				for(my $k = $#{$matrix}; $k >= 0; $k--){
					$matrix->[$k][0] =~ tr/ /_/;
					if(exists($tTpTd->{$matrix->[$k][0]}->{$infile[$f]})){
						die("Terminals must have unique names ($matrix->[$k][0] is duplicated in $infile[$f])!\n");			
						}
					for(my $j = 1; $j <= $#{$matrix->[$k]}; $j++){
						$tTpTd->{$matrix->[$k][0]}->{$infile[$f]} .= $matrix->[$k][$j];
						}
					}
				for(my $h = $fieldCounter-1; $h >= 0; $h--){
					if($add->[$f][$h] != 0){
						$add->[$f][$h] = 1;
						}
					}	
				$morpho->[$f] = 1;
				next;
				} else { ############################### READ FASTA
					if($codeIndel && !$quoteIndel){
						$quoteIndel = "\nindel characters coded using 2xread using the \"simple gap coding\" method of SIMMONS AND OCHOTERENA. 2000. Gaps as characters in sequence-based phylogenetic analysis. Systematic Biology 49: 369-381";
						}	
					my $seq = ();
					@bases = ();	
					while(my $line = <INFILE>){
						$line =~ s/[\r\n]+$//;
						if(length($line)){
							if($line =~ m/^>/){
								$line =~ tr/ \/\-/\_\_\_/;
								$line =~ tr/[A-z][0-9]\_\.//cd;
								$matrix->[$terms][0] = $line;
								if(length($seq)){
									$seq = uc($seq);
									$seq =~ tr/\?/\-/;
									$seq =~ tr/ABCDEFGHIKLMNOPQRSTUVWXYZ\-//cd;
									if($seq =~ m/E|F|I|L|O|P|Q|X|Z/){
										$aa->[$f] = 1;
										}
									@bases = split(//, $seq);
									for(my $i = $#bases; $i >= 0; $i--){
										$matrix->[$terms-1][$i+1] = $bases[$i];
										}
									}
								$seq = ();
								$terms++;	
								} else {
									$seq .= $line;
									}
							}
						}
					$seq = uc($seq);
					$seq =~ tr/\?/\-/;
					$seq =~ tr/ABCDEFGHIKLMNOPQRSTUVWXYZ\-//cd;
					if($seq =~ m/E|F|I|L|O|P|Q|X|Z/){
						$aa->[$f] = 1;
						}
					@bases = split(//, $seq);
					for(my $i = $#bases; $i >= 0; $i--){
						$matrix->[$terms-1][$i+1] = $bases[$i];
						}				
					}
		close(INFILE);
		

		
		############################### CHAR LABELS
		my $chars = $#bases + 1;
		undef(@bases);
		my @names = ();
		for(my $j = ($chars-1); $j >= 0; $j--){
			$names[$j] = $stem{$infile[$f]} . '_' .  $j;
			}
		$debug .= "$infile[$f] without indels: $terms terminals and $chars characters\n";	
		
	
		
		if($codeIndel == 1){ ############################### FIND INDELS
			my %indels = ();
			for(my $k = $#{$matrix}; $k >= 0; $k--){		
				my $lastFive = 0;
				my $lastThree = 0;
				my $inGap = 0;
				for(my $j = $chars; $j >= 1; $j--){
					if(($matrix->[$k][$j] eq '-') && ($inGap == 0)){ ### gap 3'
						$lastThree = $j;
						$inGap = 1;
						}
					if(($matrix->[$k][$j] ne '-') && ($inGap == 1)){ ### gap 5'
						$inGap = 0;
						$lastFive = $j + 1;
						if(($lastFive != 1) && ($lastThree != $chars)){ ### gaps not at ends
							$indels{$lastFive . '-' . $lastThree} = $lastFive;
							}
						}
					}
				}
			my @indelNames = sort({$indels{$a} <=> $indels{$b}} keys(%indels));	
			undef(%indels);
			for(my $j = $#indelNames; $j >= 0; $j--){
				$names[$j+$chars] = $stem{$infile[$f]} . '_indel_' . $indelNames[$j];
				}
			$debug .= 'Total number of indels: ' . ($#indelNames + 1) . "\n";	

			
			
			############################### GET ENDS
			my $ends->[0][0] = ();
			for(my $k = $#{$matrix}; $k >= 0; $k--){
				for(my $j = 1; $j <= $chars; $j++){
					if($matrix->[$k][$j] =~ m/A|B|C|D|E|F|G|H|I|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z/){
						$ends->[$k][0] = $j;
						last;
						}
					}	
				for(my $j = $chars; $j >= 1; $j--){
					if($matrix->[$k][$j] =~ m/A|B|C|D|E|F|G|H|I|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z/){
						$ends->[$k][1] = $j;
						last;
						}	
					}			
				}
			
			
				
			############################### SCORE INDELS	
			for(my $i = $#indelNames; $i >= 0; $i--){
				(my $fiveEnd, my $threeEnd) = split(/-/, $indelNames[$i]);
				for(my $k = $#{$matrix}; $k >= 0; $k--){		
					my $score = '?';
					if(($fiveEnd > $ends->[$k][0]) && ($threeEnd < $ends->[$k][1])){
						for(my $j = $fiveEnd; $j <= $threeEnd; $j++){
							if($matrix->[$k][$j] =~ m/A|B|C|D|E|F|G|H|I|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z/){
								$score = 0;
								last;
								}
							}
						if($score eq '?'){
							if(($matrix->[$k][$fiveEnd-1] ne '-') && ($matrix->[$k][$threeEnd+1] ne '-')){
								$score = 1;
								} else {
									$score = '-';
									}
							}
						}
					$matrix->[$k][$chars + $i + 1] = $score;	
					}
				}	
				$partIndels->[$f] = $#indelNames + 1;
			}
		$partSize->[$f] = $chars;

		
		############################### FORMAT DATA
		my $charsCounter = $#{$matrix->[$#{$matrix}]};
		for(my $k = $#{$matrix}; $k >= 0; $k--){
			if(exists($tTpTd->{$matrix->[$k][0]}->{$infile[$f]})){
				die("Terminals must have unique names ($matrix->[$k][0] is duplicated in $infile[$f])!\n");			
				}
			if($charsCounter != $#{$matrix->[$k]}){
				die("Terminals do not have the same number of characters in matrix $infile[$f]!\n");
				}
			for(my $j = 1; $j <= $#{$matrix->[$k]}; $j++){
				$tTpTd->{$matrix->[$k][0]}->{$infile[$f]} .= $matrix->[$k][$j];
				}
			}
		@{$pTc->{$infile[$f]}} = @names;
		$debug .= "$infile[$f] with indels: $terms terminals and $partSize->[$f] + $partIndels->[$f] characters\n";
		}
	@taxa = sort({$b cmp $a} keys(%{$tTpTd}));
	$debug .= "Number of taxa: $#taxa\n";
	my $max = 0;
	for(my $k = $#taxa; $k >= 0; $k--){					
		if($max < length($taxa[$k])){
			$max = length($taxa[$k]);
			}
		}
	$max += 3;

	$quote .= $quoteIndel;
	
	###############################	PRINT PRETTY
	my $buffer = ();
	if($xread == 1){
		############################### xread preample
		for(my $r = $#infile; $r >= 0; $r--){
			if($aa->[$r]){
				$buffer = "nstates 32;\n";
				last;
				}
			}
		$buffer .= "xread\n'$quote'\n";
		my $chars = ();
		for(my $r = $#infile; $r >= 0; $r--){
			$chars += $partIndels->[$r];
			$chars += $partSize->[$r];
			}
		$buffer .= $chars . ' ' . ($#taxa + 1) . "\n";
			
		############################### print data
		open(FILE,'>',"$outputName.ss") or die("Cannot write to $outputName.ss!\n");
		for(my $q = $#taxa; $q >= 0; $q--){			
			$buffer .= $taxa[$q] . (' ' x ($max - length($taxa[$q])));
			for(my $r = $#infile; $r >= 0; $r--){
				if(!exists($tTpTd->{$taxa[$q]}->{$infile[$r]})){
					$buffer .= '?' x ($partIndels->[$r] + $partSize->[$r]);
					} else {
						my $residue = $tTpTd->{$taxa[$q]}->{$infile[$r]};
						if($aa->[$r]){
							$residue =~ s/B/[DN]/;
							$residue =~ s/X/[ACDEFGHIKLMNOPQRSTUVWY]/;
							$residue =~ s/Z/[EQ]/;
							$residue =~ tr/ACDEFGHIKLMNOPQRSTUVWY/0123456789ABCDEFGHIJKL/;
							} else {
								$residue =~ s/R/[AG]/g;
								$residue =~ s/Y/[CT]/g;
								$residue =~ s/M/[AC]/g;
								$residue =~ s/K/[GT]/g;
								$residue =~ s/S/[CG]/g;
								$residue =~ s/W/[AT]/g;
								$residue =~ s/H/[ACT]/g;
								$residue =~ s/B/[CGT]/g;
								$residue =~ s/V/[ACG]/g;
								$residue =~ s/D/[AGT]/g;
								$residue =~ s/N/[ACGT]/g;
								$residue =~ tr/ACGTU/01233/;
								}
						$buffer .= $residue;
						if(length($buffer > 10000)){
							print(FILE $buffer);
							$buffer = ();
							}
						}
				}
			$buffer .= "\n";
			}

		############################### chars
		$buffer .= ";\ncc - .;\n";
		my $adds = ();
		my $j = 0;
		for(my $r = $#infile; $r >= 0; $r--){
			if($morpho->[$r] == 1){
				for(my $k = 0; $k < $partSize->[$r]; $k++){
					if($add->[$r][$k] == 1){
						$adds .= ' + ' . ($k + $j);
						}
					}
				$j += $partSize->[$r];	
				} else {
					$j += $partSize->[$r] + $partIndels->[$r];
					}
			}
		if($adds){
			$buffer .= 'cc' . $adds . ";\n";
			}
		
		$buffer .= "proc/;\n\n#\n\$\n;\ncn ";
		$j = 0;
		for(my $r = $#infile; $r >= 0; $r--){
			for(my $k = 0; $k <= $#{$pTc->{$infile[$r]}}; $k++){
				$buffer .= '{' . $j . ' ' . ${$pTc->{$infile[$r]}}[$k];
				if(${$pTc->{$infile[$r]}}[$k] =~ m/.* .* .*/){
					$buffer .= "\n";
					} elsif((${$pTc->{$infile[$r]}}[$k] =~ m/indel/)){
						$buffer .= " absent present /;\n";
						} else {
							if($aa->[$r]){
								$buffer .= " A C D E F G H I K L M N O P Q R S T U V W Y /;\n";
								} else {
									$buffer .= " A C G T /;\n";
									}
							}
				$j++;		
				}
			}

		print(FILE $buffer);
		}
		
		
		
	if($phylip == 1){ ############################### phylip
		############################### create partition file
		my $chars = 0;
		my $partitions = ();
		my $partCounter = 1;
		my $additiveWarning = 0;
		for(my $r = $#infile; $r >= 0; $r--){
			if($morpho->[$r] == 1){
				$partitions .= 'MULTI, p' . $partCounter . '=' . ($chars + 1) . '-' . ($chars + $partSize->[$r]) . "\n";
				$chars += $partSize->[$r];
				$partCounter += 1;
				if($additiveWarning == 0){
					for(my $c = $#{$add->[$r]}; $c >= 0; $c--){
						if($add->[$r][$c] == 1){
							print(STDERR "Additive characters have been changed to non-additive in extended PHYLIP files (RAxML does not support additivity).\n");
							$additiveWarning = 1;
							last;
							}
						}
					}
				} elsif(($codeIndel == 1) && ($partIndels->[$r] >= 1)){
					if($aa->[$r]){
						$partitions .= 'WAG, p';
						} else {
							$partitions .= 'DNA, p';
							}
					$partitions .= $partCounter . '=' . ($chars + 1) . '-' . ($chars + $partSize->[$r]) . "\nBIN, p" . ($partCounter + 1) . '=' . ($chars + $partSize->[$r] + 1) . '-' . ($chars + $partSize->[$r] + $partIndels->[$r]) . "\n";
					$chars += $partIndels->[$r];
					$chars += $partSize->[$r];
					$partCounter += 2;
					} else {
						if($aa->[$r]){
							$partitions .= 'WAG, p';
							} else {
								$partitions .= 'DNA, p';
								}
						$partitions .= $partCounter . '=' . ($chars + 1) . '-' . ($chars + $partSize->[$r]) . "\n";
						$chars += $partSize->[$r];
						$partCounter += 1;
						}
			}
		open(FILE,'>', "$outputName.part");
		print(FILE $partitions);
		close(FILE);
	
		############################### print data
		open(FILE,'>',"$outputName.phy") or die("Cannot write to $outputName.phy!\n");
		$buffer = ($#taxa + 1) . ' ' . $chars . "\n";	
		for(my $q = $#taxa; $q >= 0; $q--){
			$buffer .= $taxa[$q] . (' ' x ($max - length($taxa[$q])));
			for(my $r = $#infile; $r >= 0; $r--){
				if(!exists($tTpTd->{$taxa[$q]}->{$infile[$r]})){
					$buffer .= '?' x ($partIndels->[$r] + $partSize->[$r]);
					} elsif($tTpTd->{$taxa[$q]}->{$infile[$r]} =~ m/\[/){
						my $things = $tTpTd->{$taxa[$q]}->{$infile[$r]};
						$things  =~ s/\[\d+\]/?/g;
						if($RAXMLpoly == 0){
							print(STDERR "Polymorphisms have been changed to missing data in RAxML configuration file (RAxML does not accept polymorphisms for multistate characters).\n");
							$RAXMLpoly = 1;
							}
						$buffer .= $things;
						} else {
							$buffer .= $tTpTd->{$taxa[$q]}->{$infile[$r]};
							}
				}
			$buffer .= "\n";
			if(length($buffer > 10000)){
				print(FILE $buffer);
				$buffer = ();
				}
			}
		print(FILE $buffer);
		close(FILE);
		}

		
		
	if($nexus == 1){
		############################### preample
		my $partitions = ();
		my $chars = 0;
		my $terms = ($#taxa + 1);
		my $garliBuffer = "#nexus\n[$quote]\nbegin taxa;\n\tdimensions ntax=$terms;\n\ttaxlabels\n";
		my $mrBayesBuffer = "#nexus\n[$quote]\nbegin taxa;\n\tdimensions ntax=$terms;\n\ttaxlabels\n";
		my $setsBlock = ();
		my $dataType = '(';
		my $ordModels = ();

		###############################  sets block
		$setsBlock .= "begin sets;\n";
		my $chunks = ();
		my $partCounter = 1;
		for(my $r = $#infile; $r >= 0; $r--){
			if($morpho->[$r] == 1){
				$setsBlock .= "\tcharset Part" . $partCounter . ' = ' . ($chars + 1) . '-' . ($chars + $partSize->[$r]) . "\n";
				$dataType .= 'standard:' . ($chars + 1) . '-' . ($chars + $partSize->[$r]);
				$chars += $partSize->[$r];
				$chunks .= 'chunk' . $partCounter . ': Part' . $partCounter;
				$partCounter++;
				if($r > 0){
					$chunks .= ', ';
					$dataType .= ','
					} else{
						$chunks .= ';';
						}
				} elsif(($codeIndel == 1) && ($partIndels->[$r] >= 1)){
					$setsBlock .= "\tcharset Part" . $partCounter . ' = ' . ($chars + 1) . '-' . ($chars + $partSize->[$r]) . "\n\tcharset Part" . ($partCounter + 1) . ' = ' . ($chars + $partSize->[$r] + 1) . '-' . ($chars + $partSize->[$r] + $partIndels->[$r]) . "\n";
					if($aa->[$r]){
						$dataType .= 'protein:';
						}else{
							$dataType .= 'DNA:';
							}
					$dataType .= ($chars + 1) . '-' . ($chars + $partSize->[$r]) . ',standard:' . ($chars + $partSize->[$r] + 1) . '-' . ($chars + $partSize->[$r] + $partIndels->[$r]);
					$chars += ($partIndels->[$r] + $partSize->[$r]);
					$chunks .= 'chunk' . $partCounter . ': Part' . $partCounter . ', chunk' . ($partCounter + 1) . ': Part' . ($partCounter + 1);
					$partCounter += 2;
					if($r > 0){
						$chunks .= ', ';
						$dataType .= ','
						} else{
							$chunks .= ';';
							}
					} else{
						$setsBlock .= "\tcharset Part" . $partCounter . ' = ' . ($chars + 1) . '-' . ($chars + $partSize->[$r]) . "\n";
						if($aa->[$r]){
							$dataType .= 'protein:';
							}else{
								$dataType .= 'DNA:';
								}
						$dataType .= ($chars + 1) . '-' . ($chars + $partSize->[$r]);
						$chars += $partSize->[$r];
						$chunks .= 'chunk' . $partCounter . ': Part' . $partCounter;
						$partCounter++;
						if($r > 0){
							$chunks .= ', ';
							$dataType .= ','
							} else{
								$chunks .= ';';
								}
						}
			}
		$setsBlock .= "\tcharpartition default = $chunks\nend;\n";
		$dataType .= ')';

		############################### get max label size
		$max = 0;
		for(my $k = $#taxa; $k >= 0; $k--){
			$garliBuffer .= "\t\t$taxa[$k]\n";
			$mrBayesBuffer .= "\t\t$taxa[$k]\n";
			if($max < length($taxa[$k])){
				$max = length($taxa[$k]);
				}
			}
		$garliBuffer .= "\t\t;\nend;\n";	
		$mrBayesBuffer .= "\t\t;\nend;\n";	
		$max += 3;

		my $charsBlock = ();
		my $charsBlockAdditive = ();
		my %symbols = ();
		my %symbolsAdditive = ();
		my $jj = 0;
		for(my $r = $#infile; $r >= 0; $r--){
			my $j = 0;
			my $i = 0;
			if($morpho->[$r] == 1){ ############################### get non-molecular data from partition $r
				#### insert character names
				my $charLabels = "\tcharlabels\n";
				my $stateLabels = "\tstatelabels\n";
				my $charLabelsAdditive = "\tcharlabels\n";
				my $stateLabelsAdditive = "\tstatelabels\n";
				for(my $k = 0; $k <= $#{$pTc->{$infile[$r]}}; $k++){
					${$pTc->{$infile[$r]}}[$k] =~ tr/;//d;
					my @bits = split(/ /,${$pTc->{$infile[$r]}}[$k]);
					if($add->[$r][$k] == 0){
						$j++;
						$jj++;
						$charLabels .= "\t\t[$j] $bits[0] \n";
						$stateLabels .= "\t\t$j";
						for(my $l = 1; $l <= $#bits; $l++){
							$stateLabels .= " $bits[$l]";
							}
						$stateLabels .= ",\n";
						} elsif($add->[$r][$k] == 1){
							$i++;
							$jj++;
							$charLabelsAdditive .= "\t\t[$i] $bits[0] \n";
							$stateLabelsAdditive .= "\t\t$i";
							for(my $l = 1; $l <= $#bits; $l++){
								$stateLabelsAdditive .= " $bits[$l]";
								}
							$stateLabelsAdditive .= ",\n";
							}
					}
				$charsBlock = "begin characters;\n\tdimensions nchar=$j;\n\tformat datatype=standard symbols=\"this-is-where-the-list-of-symbols-will-added-by-2matrix\" missing=? gap=-;\n" . $charLabels . "\t\t;\n" . $stateLabels . "\t\t;\n\tmatrix\n";
				$charsBlockAdditive = "begin characters;\n\tdimensions nchar=$i;\n\tformat datatype=standard symbols=\"this-is-where-the-list-of-symbols-will-added-by-2matrix\" missing=? gap=-;\n" . $charLabelsAdditive . "\t\t;\n" . $stateLabelsAdditive . "\t\t;\n\tmatrix\n";
				} else{ #### insert character names
					my $charLabels = "\tcharlabels\n";
					my $stateLabels = "\tstatelabels\n";
					for(my $k = 0; $k < $partSize->[$r]; $k++){
						${$pTc->{$infile[$r]}}[$k] =~ tr/;//d;
						my @bits = split(/ /,${$pTc->{$infile[$r]}}[$k]);
						if($add->[$r][$k] == 0){
							$j++;
							$jj++;
							$charLabels .= "\t\t[$j] $bits[0] \n";
							if($aa->[$r]){
								$stateLabels .= "\t\t$j A C D E F G H I K L M N O P Q R S T U V W Y,\n";
								} else {
									$stateLabels .= "\t\t$j A C G T,\n";
									}
							}
						}
					if($aa->[$r]){
						$charsBlock = "begin characters;\n\tdimensions nchar=$partSize->[$r];\n\tformat datatype=protein missing=? gap=-;\n" . $charLabels . "\t\t;\n" . $stateLabels . "\t\t;\n\tmatrix\n";
						} else {
							$charsBlock = "begin characters;\n\tdimensions nchar=$partSize->[$r];\n\tformat datatype=DNA missing=? gap=-;\n" . $charLabels . "\t\t;\n" . $stateLabels . "\t\t;\n\tmatrix\n";
							}
					$j = $partSize->[$r];
					}
			for(my $q = $#taxa; $q >= 0; $q--){
				$charsBlock .= "\t\t" . $taxa[$q];
				$charsBlock .= ' ' x ($max - length($taxa[$q]));
				$charsBlockAdditive .= "\t\t" . $taxa[$q];
				$charsBlockAdditive .= ' ' x ($max - length($taxa[$q]));
				if(!exists($tTpTd->{$taxa[$q]}->{$infile[$r]})){
					$charsBlock .= '?' x $j;
					$charsBlockAdditive .= '?' x $i;
					} elsif($morpho->[$r] == 1){
						my @data = split(//, $tTpTd->{$taxa[$q]}->{$infile[$r]});
						my $c = 0;
						for(my $k = 0; $k <= $#data; $k++){
							my $state = ();
							if($data[$k] eq '['){
								while($data[$k] ne ']'){
									$state .= $data[$k];
									if($data[$k] =~ m/\d/){
										if($add->[$r][$c] == 1){
											$symbolsAdditive{$data[$k]} = 1;
											} else {
												$symbols{$data[$k]} = 1;
												}
										}
									$k++;
									}
								$state .= ']';	
								} else {
									$state = $data[$k];
									if($data[$k] =~ m/\d/){
										if($add->[$r][$c] == 1){
											$symbolsAdditive{$state} = 1;
											} else {
												$symbols{$state} = 1;
												}
										}
									}
							$state =~ tr/[]/{}/;
							$state =~ tr/[0-9]{},\?\-//cd;
							if($add->[$r][$c] == 1){
								$charsBlockAdditive .= $state;
								} else {
									$charsBlock .= $state;
									}
							$c++;	
							}	
						} else { #$part = substr($string, $initial_position, $length);
							$charsBlock .= substr($tTpTd->{$taxa[$q]}->{$infile[$r]}, 0, $partSize->[$r]);
							}
				$charsBlock .= "\n";
				$charsBlockAdditive .= "\n";
				}
			$charsBlock .= "\t\t;\nend;\n";
			$charsBlockAdditive .= "\t\t;\nend;\n";
			my $add = join('', sort({$a <=> $b} keys(%symbols)));
			$charsBlock =~ s/symbols="this-is-where-the-list-of-symbols-will-added-by-2matrix"/symbols="$add"/;
			$add = join('', sort({$a <=> $b} keys(%symbolsAdditive)));
			$charsBlockAdditive =~ s/symbols="this-is-where-the-list-of-symbols-will-added-by-2matrix"/symbols="$add"/;
			if($j == $partSize->[$r]){
				$garliBuffer .= $charsBlock;
				$ordModels->[$r] = 0;
				} elsif($i == $partSize->[$r]){
					$garliBuffer .= $charsBlockAdditive;
					$ordModels->[$r] = 1;
					} else {
						$garliBuffer .= $charsBlock . "\n" . $charsBlockAdditive;
						$ordModels->[$r] = 2;
						}
			undef($charsBlockAdditive);	
			$charsBlock = ();
			if(($codeIndel == 1) && ($partIndels->[$r] >= 1)){
				$j = 0;
				#### insert character names
				my $charLabels = "\tcharlabels\n";
				my $stateLabels = "\tstatelabels\n";
				for(my $k = ($partSize->[$r]); $k < ($partSize->[$r] + $partIndels->[$r]); $k++){
					${$pTc->{$infile[$r]}}[$k] =~ tr/;//d;
					${$pTc->{$infile[$r]}}[$k] =~ s/\-/_to_/g;
					my @bits = split(/ /,${$pTc->{$infile[$r]}}[$k]);
					$j++;
					$jj++;
					$charLabels .= "\t\t[$j] $bits[0] \n";
					$stateLabels .= "\t\t$j absent present,\n";
					}

				############################### get indel data from partition $r
				$charsBlock = "begin characters;\n\tdimensions nchar=$partIndels->[$r];\n\tformat datatype=standard symbols=\"01\" missing=? gap=-;\n" . $charLabels . "\t\t;\n" . $stateLabels . "\t\t;\n\tmatrix\n";
				for(my $q = $#taxa; $q >= 0; $q--){
					$charsBlock .= "\t\t" . $taxa[$q];
					$charsBlock .= ' ' x ($max - length($taxa[$q]));
				if(!exists($tTpTd->{$taxa[$q]}->{$infile[$r]})){
					for(my $z = $partIndels->[$r]; $z >= 1; $z--){
						$charsBlock .= '?';
						}
					} else{ ### $part = substr($string, $initial_position, $length);
						$charsBlock .= substr($tTpTd->{$taxa[$q]}->{$infile[$r]}, $partSize->[$r], $partIndels->[$r]);
						}
					$charsBlock .= "\n";
					}
				$charsBlock .= "\t\t;\nend;\n";
				$garliBuffer .= $charsBlock;
				$charsBlock = ();
				}
			}
		my $dataBlock = "begin data;\n\tdimensions ntax=" . ($#taxa + 1) . " nchar=$chars;\n\tformat datatype=mixed$dataType missing=? gap=-;\n\tmatrix\n";
		for(my $q = $#taxa; $q >= 0; $q--){
			for(my $r = $#infile; $r >= 0; $r--){
				############################### get nucleotide data from partition $r
				if($r == $#infile){
					$dataBlock .= "\t\t" . $taxa[$q];
					$dataBlock .= ' ' x ($max - length($taxa[$q]));
					}
				if(!exists($tTpTd->{$taxa[$q]}->{$infile[$r]})){
					for(my $z = ($partSize->[$r] + $partIndels->[$r]); $z >= 1; $z--){
						$dataBlock .= '?';
						}
					} else{
						$dataBlock .= $tTpTd->{$taxa[$q]}->{$infile[$r]};
						}
				}
			$dataBlock .= "\n";
			}
		$dataBlock =~ tr/[]/{}/;
		$dataBlock .= "\t\t;\nend;\n";
		$mrBayesBuffer .= $dataBlock;
		$dataBlock = ();

		$mrBayesBuffer .= $setsBlock;

		my $adds = ();
		my $j = 0;
		for(my $r = $#infile; $r >= 0; $r--){
			if($morpho->[$r] == 1){
				for(my $k = 0; $k < $partSize->[$r]; $k++){
					if($add->[$r][$k] == 1){
						$adds .= ($k + $j + 1) . ' ';
						}
					}
				$j += $partSize->[$r];	
				} else {
					$j += $partSize->[$r] + $partIndels->[$r];
					}
			}
		
		my $mrBayesBlock = "\nbegin mrbayes;\n	log start filename= $outputName.log;\n";
		if($adds){
			$mrBayesBlock .= "\tctype ordered: $adds;\n";
			}
		$mrBayesBlock .= "\tset autoclose=yes;\n\tlog stop;\nend;\n\n";
		$mrBayesBuffer .= $mrBayesBlock;
	
		open(FILE,'>',"$outputName.garli.nex") or die("Could not open $outputName.garli.nex\n");
		print(FILE $garliBuffer);	
		close(FILE);	
		$garliBuffer = ();
		
		open(FILE,'>',"$outputName.mrbayes.nex") or die("Could not open $outputName.mrbayes.nex\n");
		print(FILE $mrBayesBuffer);	
		close(FILE);	
		$mrBayesBuffer = ();
		
		$garliBuffer = "[general]\ndatafname = $outputName.garli.nex\nconstraintfile = none\nstreefname = random\nattachmentspertaxon = 100\nofprefix = mixedDnaMkv\nrandseed = -1\navailablememory = 512\nlogevery = 10\nsaveevery = 100\nrefinestart = 1\noutputeachbettertopology = 0\noutputcurrentbesttopology = 0\nenforcetermconditions = 1\ngenthreshfortopoterm = 10000\nscorethreshforterm = 0.001\nsignificanttopochange = 0.01\noutputphyliptree = 0\noutputmostlyuselessfiles = 0\nwritecheckpoints = 0\nrestart = 0\noutgroup = 1\nresampleproportion = 1.0\ninferinternalstateprobs = 0\noutputsitelikelihoods = 0\noptimizeinputonly = 0\ncollapsebranches = 1\nsearchreps = 5\nbootstrapreps = 0\nlinkmodels = 0\nsubsetspecificrates = 1\n\n";
		$partCounter = 0;
		for(my $r = $#infile; $r >= 0; $r--){
			if($morpho->[$r] == 1){
				if($ordModels->[$r] == 0){
					$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = standard\nratematrix = 1rate\nstatefrequencies = equal\nratehetmodel = none\nnumratecats = 1\ninvariantsites = none\n\n";
					$partCounter++;
					} elsif($ordModels->[$r] == 1){
						$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = standardordered\nratematrix = 1rate\nstatefrequencies = equal\nratehetmodel = none\nnumratecats = 1\ninvariantsites = none\n\n";
						$partCounter++;
						} elsif($ordModels->[$r] == 2) { 
							$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = standard\nratematrix = 1rate\nstatefrequencies = equal\nratehetmodel = none\nnumratecats = 1\ninvariantsites = none\n\n[model" . ($partCounter + 2) . "]\ndatatype = standardordered\nratematrix = 1rate\nstatefrequencies = equal\nratehetmodel = none\nnumratecats = 1\ninvariantsites = none\n\n";
							$partCounter += 2;
							}
				} elsif(($codeIndel == 1) && ($partIndels->[$r] >= 1)){
					if($aa->[$r]){
						$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = aminoacid\nratematrix = WAG\n";
						}else{
							$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = nucleotide\nratematrix = 6rate\n";
							}
					$garliBuffer .= "statefrequencies = estimate\nratehetmodel = gamma\nnumratecats = 4\ninvariantsites = estimate\n\n[model" . ($partCounter + 2) . "]\ndatatype = standard\nratematrix = 1rate\nstatefrequencies = equal\nratehetmodel = none\nnumratecats = 1\ninvariantsites = none\n\n";
					$partCounter += 2;
					} else{
						if($aa->[$r]){
							$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = aminoacid\nratematrix = WAG\n";
							}else{
								$garliBuffer .= '[model'. ($partCounter + 1) . "]\ndatatype = nucleotide\nratematrix = 6rate\n";
								}
						$garliBuffer .= "statefrequencies = estimate\nratehetmodel = gamma\nnumratecats = 4\ninvariantsites = estimate\n\n";
						$partCounter++;
						}
			}
		$garliBuffer .= "[master]\nnindivs = 4\nholdover = 1\nselectionintensity = 0.5\nholdoverpenalty = 0\nstopgen = 5000000\nstoptime = 5000000\nstartoptprec = 0.5\nminoptprec = 0.01\nnumberofprecreductions = 10\ntreerejectionthreshold = 50.0\ntopoweight = 0.01\nmodweight = 0.002\nbrlenweight = 0.002\nrandnniweight = 0.1\nrandsprweight = 0.3\nlimsprweight =  0.6\nintervallength = 100\nintervalstostore = 5\nlimsprrange = 6\nmeanbrlenmuts = 5\ngammashapebrlen = 1000\ngammashapemodel = 1000\nuniqueswapbias = 0.1\ndistanceswapbias = 1.0\n\n";
		open(FILE,'>',"$outputName.conf") or die("Could not open $outputName.conf");
		print(FILE $garliBuffer);
		close(INFILE);
		$garliBuffer = ();
		}

		
		
	} else {
		############################### PRINT INFO
		print("\nA PERL script for merging and translating FASTA alignments into XREAD, extended\n");
		print("PHYLIP and NEXUS formats with indel characters coded using the 'simple' gap\n");
		print("coding method of Simmons and Ochoterena (2000; Gaps as characters\n");
		print("in sequence-based phylogenetic analysis. Systematic Biology 49: 369-381).\n\n");
		print("USAGE: 2matrix.pl -i <infile_1> [ -i <infile_2>... ] -n <root-name> [ -d ]\n");
		print("\t -o n|x|p [ -s <stem-name_1> [ -s <stem-name_2> ] ... ]\n\n");
		print("WHERE:\n\n\t-d\tDo NOT code indels.\n\n");
		print("\t-i\tSpecifies a matrix (aligned FASTA, csv, or xread cf.\n");
		print("\t\tHennig86/NONA/WinClada). If several matrices are to be merged,\n");
		print("\t\tfilenames should be input with multiple -i flags.\n\n");
		print("\t-n\t<root-name> for output files.\n\n");
		print("\t-o\tSets the output format: 'x' for XREAD, 'n' for NEXUS,\n\t\tand 'p' for extended PHYLIP. If PHYLIP format is selected,\n");
		print("\t\ta RAxML partition file will automatically be created\n");
		print("\t\t(<root-name>.part). If NEXUS format is selected, files\n");
		print("\t\tcompatible with both Garli (<root-name>.garli.nex) and MrBayes\n");
		print("\t\t(<root-name>.mrbayes.nex) will be created. Additionally, a\n");
		print("\t\tGarli configuration file will be automatically generated\n\t\t(<root-name>.conf).\n\n");
		print("\t-s\t<stem-name> for naming sequence characters (XREAD and NEXUS\n\t\tonly).\n\n"); 
		print("If you use this program please cite Salinas and Little (2014; 2matrix: a utility\n");
		print("for indel coding and phylogenetic matrix concatenation. Applications in Plant\nSciences 2: 1300083).\n\n");
		}
		
		
		
exit(0);
