package bio.uima.util;

import java.util.Map;

import com.json.parsers.JSONParser;
import com.json.parsers.JsonParserFactory;

public abstract class MatchScoring {

	public static final byte LINEAR_PENALTY 	= 0;
	public static final byte CONSTANT_PENALTY 	= 1;
	public static final byte AFFINE_PENALTY 	= 2;
	public static final byte CONVEX_PENALTY 	= 3;
	public static final byte PROFILE_PENALTY 	= 4;
	
	private byte gapType;
	private double gapConst;
	private double gapLinear;
	private IGapProfile profile;
	
	private Map subMatrix;
	private double subsConstant;
	private boolean isLikelihood;
	
	public MatchScoring(byte gapPenalty) {
		switch(gapType) {
			case CONSTANT_PENALTY:
			case LINEAR_PENALTY:
			case AFFINE_PENALTY:
			case CONVEX_PENALTY:
			case PROFILE_PENALTY:
				break;
			default:
				throw new IllegalArgumentException("Unknown gap penalty type: " + gapType);
		}
		gapType = gapPenalty;
	}
	
	public void setGapConstant(double c) {
		gapConst = c;
	}
	
	public double getGapConstant() {
		return gapConst;
	}
	
	public void setGapProfile(IGapProfile p) {
		profile = p;
	}
	
	public IGapProfile getGapProfile() {
		return profile;
	}
	
	public void setGapLinear(double l) {
		gapLinear = l;
	}
	
	public double getGapLinear() {
		return gapLinear;
	}
	
	public double gapPenalty(int length) {
		switch(gapType) {
			case CONSTANT_PENALTY:
				return gapConst;
			case LINEAR_PENALTY:
				return gapLinear*length;
			case AFFINE_PENALTY:
				return gapConst + gapLinear*length;
			case CONVEX_PENALTY:
				return gapConst + gapLinear*Math.log(length);
			case PROFILE_PENALTY:
				if (profile == null) throw new IllegalStateException("Cannot calculate profile penalty.  Profile has not been initialized.");
				return profile.getPenalty(length);
			default:
				throw new IllegalArgumentException("Unknown gap penalty type: " + gapType);
		}
	}
	
	public void setSubstitutionMatrix(Map m) {
		subMatrix = m;
	}
	
	public void setSubstitutionMatrix(String json) {
		JsonParserFactory factory=JsonParserFactory.getInstance();
		JSONParser parser=factory.newJsonParser();
		subMatrix=parser.parseJson(json);
		isLikelihood = "likelihood".equals((String)subMatrix.get("matrix-type"));
	}
	
	public Map getSubstitutionMatrix() {
		return subMatrix;
	}
	
	public void setSubsConstant(double d) {
		subsConstant = d;
	}
	public double getSubsConstant() {
		return subsConstant;
	}
	public void setSubMatrixIsLikelihood(boolean b) {
		isLikelihood = b;
	}
	
	public double substituteCost(String c1, String c2) {
		if (subMatrix == null) return subsConstant;
		try { 
			return (isLikelihood ? -1 : 1)*Double.parseDouble((String)((Map)subMatrix.get(c1)).get(c2));
		}
		catch (NullPointerException ex) {
			NullPointerException nex = new NullPointerException("Could not find substitution matrix entry for " + c1 + "->" +c2);
			nex.setStackTrace(ex.getStackTrace());
			throw nex;
		}
	}
}
