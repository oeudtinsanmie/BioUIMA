package bio.uima;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;
import org.apache.uima.resource.ResourceAccessException;

import bio.uima.util.DNAReader;

import com.json.generators.JSONGenerator;
import com.json.generators.JsonGeneratorFactory;

/*
 * This annotator should read DNA sequences from the CAS
 * and translate each sequence to protein sequences.  
 * 
 * See http://en.wikipedia.org/wiki/DNA_codon_table
 *     http://en.wikipedia.org/wiki/Open_reading_frame
 */
public class ProteinSequenceAnnotator extends JCasAnnotator_ImplBase  {
	private String[] proteinWords = new String[6];
	
	private char[][] frames = new char[proteinWords.length][3];
	private int head=0;
	private Vector proteins;
	
	private Map data = new HashMap();
	private int[] transcribing = new int[proteinWords.length];
	
	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		// read DNA sequences from the CAS 
//		String[] seqs = cas.getDocumentText().split(" ");
//		System.out.println(seqs[0]);
//		System.out.println(seqs[1]);
		Vector seqs = new Vector();
		
		try {
			InputStream stream = getContext().getResourceAsStream("dna");
			BufferedReader buf = new BufferedReader(new InputStreamReader(stream));
			String line; int comment;
			try {
				line = buf.readLine();
				while (line != null) {
					comment = line.indexOf('#');
					if (comment >=0) {
						line = line.substring(comment);
					}
					if (!"".equals(line)) {
						seqs.add(line);
					}

					line = buf.readLine();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		} catch (ResourceAccessException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		// TODO: implement this
		// translate each DNA sequence to a protein sequence considering
		// all possible open reading frames and store in a CAS view
		
		Map datum;
		for (int j=0; j<seqs.size(); j++) {
			datum = new HashMap();
			proteins = new Vector();
			for (int i=0; i < proteinWords.length; i++) {
				proteinWords[i] = "";
				transcribing[i] = -1;
			}
			readDNASeq((String)seqs.get(j));
			datum.put("proteinWords", proteinWords);
			datum.put("proteins", proteins.toArray());
			data.put("sequence"+j, datum);
		} 

		JsonGeneratorFactory factory=JsonGeneratorFactory.getInstance();
        JSONGenerator generator=factory.newJsonGenerator();
		try {
			JCas orf1 = cas.createView("orf1");
			orf1.setDocumentText(generator.generateJson(data));
		} catch (CASException e) {
			e.printStackTrace();
		}			
	}

	// This method should generate the protein sequence for the specified
	// DNA sequence.  You man use regex.
	//
	// A codon is just a sequence of three characters.  Each codon maps to 
	// a protein.  For example the DNA codon 'TGT' maps to 'C' (cysteine).
	// So the sequence 'TGTTGTTGTTGT' should translate to 'CCCC'.
	//
	// See: http://en.wikipedia.org/wiki/DNA_codon_table
	private void readDNASeq(String seq) {
		frames[0][0] = seq.charAt(0);
		frames[0][1] = seq.charAt(1);
		frames[1][0] = seq.charAt(1);

		for (int i=2; i<seq.length(); i++) {
			for (int j=2; j>=0; j--) {
				head = (i-j)%3;
				frames[j][head] = seq.charAt(i);
				frames[j+3][head] = DNAReader.opposite(frames[j][head]);
				if (head == 2) {
					writeProtein(j);
					writeProtein(j+3);
				}
			}
		}
	}
	
	private void writeProtein(int j) {
		String protein = DNAReader.readCodon(frames[j][0], frames[j][1],frames[j][2]);
		if (transcribing[j] >= 0) {
			proteinWords[j] += protein;
			if ("X".equals(protein)) {
				proteins.add(((String)proteinWords[j]).substring(transcribing[j]));
				transcribing[j]  = -1;
			} 
		} else {
			if ("M".equals(protein)) {
				transcribing[j]  = proteinWords[j].length();
				proteinWords[j] += protein;
			} else {
				proteinWords[j] += protein.toLowerCase();
			}
		}
	}
}