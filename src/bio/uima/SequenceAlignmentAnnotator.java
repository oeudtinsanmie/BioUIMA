package bio.uima;

import java.util.HashMap;
import java.util.Map;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

import bio.uima.util.MatchScoring;

import com.json.generators.JSONGenerator;
import com.json.generators.JsonGeneratorFactory;
import com.json.parsers.JSONParser;
import com.json.parsers.JsonParserFactory;

/*
 * This annotator should compute a minimum cost alignment
 */
public class SequenceAlignmentAnnotator extends JCasAnnotator_ImplBase  {

	public static final byte SUBSTITUTE = 0;
	public static final byte INSERT = 1;
	public static final byte DELETE = -1;
	public static final byte INSDEL = 2;  // currently treating like INSERT, but could go back and report multiple paths in later versions
	
	private double[][] optArray;
	private int   [][] gapArray;
	private byte  [][] pathArray;
	
	private MatchScoring scoring;
	
	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		// read protein sequences from the CAS 
		try {
			String inputJsonString = cas.getView("orf1").getDocumentText();
			JsonParserFactory factory=JsonParserFactory.getInstance();
			JSONParser parser=factory.newJsonParser();
			Map sequences = parser.parseJson(inputJsonString);
			
			String[] seq1, seq2;
			Map alignData = new HashMap();
			for (int i=0; i<sequences.size(); i++) {
				for (int j=i+1; j<sequences.size(); j++) {
					seq1 = (String[])((Map)sequences.get("sequence"+i)).get("proteinWords");
					seq2 = (String[])((Map)sequences.get("sequence"+j)).get("proteinWords");
					
					for (int k=0; k<seq1.length; k++) {
						for (int m=0; m<seq2.length; m++) {
							alignData.put("seq"+i+"-"+j+"/wrd"+k+"-"+m, computeAlignment(seq1[k], seq2[m]));
						}
					}
				}
			}
			
			JsonGeneratorFactory genFactory=JsonGeneratorFactory.getInstance();
	        JSONGenerator generator=genFactory.newJsonGenerator();
			
			// TODO: implement this	
			JCas alignmentCas = cas.createView("alignment");
			alignmentCas.setDocumentText(generator.generateJson(alignData));
		} catch (CASException e) {
			e.printStackTrace();
		}
	}

	// This method should compute a minimum cost alignment 
	// for the specified sequences.  The cost of an insertion
	// or delete is 1 and the cost of substitution is 2.  
	// Return marked up strings where insertions/deletes are
	// represented with a '-'.  
	//
	// For example an alignment of:
	//
	// 'AAGT' and 'AGT' 
	//
	// could be 
	//
	// {'AAGT', 'A-GT'}
	//
	private String[] computeAlignment(String seq1, String seq2) {
		optArray  = new double[seq1.length()][seq2.length()];
		pathArray = new byte  [seq1.length()][seq2.length()];
		gapArray =  new int   [seq1.length()][seq2.length()];
		double delete;
		double insert;
		double substitute;
		double minVal;
		
		int i, j;
		
		for (i=0; i<seq1.length(); i++) {
			gapArray [i][0] = i;
			optArray [i][0] = scoring.gapPenalty(i);
			pathArray[i][0] = DELETE; // delete
		}
		for (j=0; j<seq2.length(); j++) {
			gapArray [0][j] = j;
			optArray [0][j] = scoring.gapPenalty(j);
			pathArray[0][j] = INSERT; // insert
		}
		
		// This method prefers fewer gaps, because constant, linear, affine and convex gap penalties are less or equal with fewest gaps along the path.  
		// If profiles exist that defy this rule, this method would need to have a way to check whether substitution is better than an insdel operation
		// Also, segregating inserts and deletes to prevent algorithm from combining an insert and a delete to get around a substitution cost.  
		// Interactions between gap penalty functions and substitution matrices seem unclear, but allowing both inserts and deletes seems to subvert the concept of the substitution matrices.
		for (j=0; j<seq2.length(); j++) {
			for (i=0; i<seq1.length(); i++) {
				if (seq1.charAt(i) == seq2.charAt(j)) {
					optArray [i][j] = optArray[i-1][j-1];
					pathArray[i][j] = SUBSTITUTE; // diagonal path, no substitution cost
				} else {
					if (gapArray[i-1][j] <= 0) {
						delete = optArray[i-1][j] + scoring.gapPenalty(-(gapArray[i-1][j] - 1));
					} else {
						delete = Double.MAX_VALUE;
					}
					if (gapArray[i][j-1] >= 0) {
						insert = optArray[i][j-1] + scoring.gapPenalty(gapArray[i][j-1] + 1);
					} else {
						insert = Double.MAX_VALUE;
					}
					substitute = optArray[i-1][j-1] + scoring.substituteCost(""+seq1.charAt(i), ""+seq2.charAt(j));
					
					minVal = Math.min(delete, insert);
					minVal = Math.min(minVal, substitute);
					optArray[i][j] = minVal;

					if (insert == minVal) {
						gapArray[i][j] = gapArray[i][j-1] + 1;
						if (delete == minVal) {
							pathArray[i][j] = INSDEL; 
						} else {
							pathArray[i][j] = INSERT;
						}
					} else if (delete == minVal) {
						pathArray[i][j] = DELETE;
						gapArray[i][j] = gapArray[i-1][j] - 1;
					} else {
						pathArray[i][j] = SUBSTITUTE;
						gapArray[i][j] = 0;
					}
				}
			}
		}
		
		i=seq1.length();
		j=seq2.length();
		String[] theAnswer = new String[2];
		theAnswer[0] = ""; theAnswer[1] = "";
		while (i>0 || j>0) {
			switch (pathArray[i][j]) {
				case SUBSTITUTE:
					theAnswer[0] = seq1.charAt(i) + theAnswer[0];
					theAnswer[1] = seq2.charAt(j) + theAnswer[1];
					i--; j--;
					break;
				case INSERT:
				case INSDEL:
					theAnswer[0] = "-"			  + theAnswer[0];
					theAnswer[1] = seq2.charAt(j) + theAnswer[1];
					j--;
					break;
				case DELETE:
					theAnswer[0] = seq1.charAt(i) + theAnswer[0];
					theAnswer[1] = "-" 			  + theAnswer[1];
					i--;
					break;
			}
		}

		return theAnswer;
	}
}