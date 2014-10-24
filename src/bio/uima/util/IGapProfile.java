package bio.uima.util;

import org.apache.uima.jcas.JCas;

public interface IGapProfile {

	public void setContext(JCas cas);
	public double getPenalty(int length);
}
