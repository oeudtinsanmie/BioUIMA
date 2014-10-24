package bio.uima.util;

public class DNAReader {

	public static String readCodon(char c1, char c2, char c3) {
		switch (c1) {
			case 'T':
			case 't':
				switch (c2) {
					case 'T':
					case 't':
						switch (c3) {
							case 'T':
							case 't': 
							case 'C': 
							case 'c':
								return "F";
							case 'A':
							case 'a': 
							case 'G': 
							case 'g':
								return "L";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'C':
					case 'c':
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g': 
								return "S";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'A': 
					case 'a': 
						switch (c3) {
							case 'T':
							case 't': 
							case 'C': 
							case 'c': 
								return "Y";
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g': 
								return "X";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'G': 
					case 'g': 
						switch (c3) {
							case 'T':
							case 't': 
							case 'C': 
							case 'c': 
								return "C";
							case 'A':
							case 'a': 
								return "X";
							case 'G': 
							case 'g': 
								return "W";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					default:
						throw new IllegalArgumentException("unknown DNA character: " + c2);
				}
			case 'C':
			case 'c': 
				switch (c2) {
					case 'T': 
					case 't': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "L";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'C':
					case 'c': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g': 
								return "P";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'A':
					case 'a': 
						switch (c3) {
							case 'T':
							case 't': 
							case 'C': 
							case 'c':  
								return "H";
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g': 
								return "Q";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'G':
					case 'g': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g': 
								return "R";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					default:
						throw new IllegalArgumentException("unknown DNA character: " + c2);
				}
			case 'A': 
			case 'a': 
				switch (c2) {
					case 'T': 
					case 't': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c':
							case 'A': 
							case 'a':
								return "I"; 
							case 'G': 
							case 'g': 
								return "M";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'C':
					case 'c': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g': 
								return "T";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'A':
					case 'a':
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c':
								return "N";
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "K";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'G':
					case 'g': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c':
								return "S";
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "R";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					default:
						throw new IllegalArgumentException("unknown DNA character: " + c2);
				}
			case 'G': 
				switch (c2) {
					case 'T':
					case 't': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "V";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'C':
					case 'c': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "A";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'A':
					case 'a': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c':
								return "D";
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "E";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					case 'G':
					case 'g': 
						switch (c3) {
							case 'T': 
							case 't': 
							case 'C': 
							case 'c': 
							case 'A': 
							case 'a': 
							case 'G': 
							case 'g':
								return "G";
							default:
								throw new IllegalArgumentException("unknown DNA character: " + c3);
						}
					default:
						throw new IllegalArgumentException("unknown DNA character: " + c2);
						
				}
			default:
				throw new IllegalArgumentException("unknown DNA character: " + c1);
		}
	}
	
	public static char opposite(char c) {
		switch(c) {
			case 'T':
			case 't':
				return 'A';
			case 'A':
			case 'a':
				return 'T';
			case 'C':
			case 'c':
				return 'G';
			case 'G':
			case 'g':
				return 'C';
			default:
				throw new IllegalArgumentException("unknown DNA character: " + c);
		}
	}
}
