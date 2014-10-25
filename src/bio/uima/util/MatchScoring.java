package bio.uima.util;

import java.util.HashMap;
import java.util.Map;

import com.json.generators.JSONGenerator;
import com.json.generators.JsonGeneratorFactory;
import com.json.parsers.JSONParser;
import com.json.parsers.JsonParserFactory;

public class MatchScoring {

	public static final byte LINEAR_PENALTY 	= 0;
	public static final byte CONSTANT_PENALTY 	= 1;
	public static final byte AFFINE_PENALTY 	= 2;
	public static final byte CONVEX_PENALTY 	= 3;
	public static final byte PROFILE_PENALTY 	= 4;
	
	private static final String LINEAR_PENALTY_STR		= "linear";
	private static final String CONSTANT_PENALTY_STR 	= "constant";
	private static final String AFFINE_PENALTY_STR 		= "affine";
	private static final String CONVEX_PENALTY_STR 		= "convex";
	private static final String PROFILE_PENALTY_STR 	= "profile";

	private static final String GAP_KEY 		= "gapType";
	private static final String GAPCONST_KEY 	= "gapConst";
	private static final String GAPLIN_KEY 		= "gapLinear";
	private static final String SUBCONST_KEY 	= "subsConst";
	private static final String SUBLIKLI_KEY 	= "subsIsLikelihood";
	private static final String SUBMAP_KEY		= "subsMap";
	
	private byte gapType;
	private double gapConst;
	private double gapLinear;
	private IGapProfile profile;
	
	private Map subMatrix;
	private double subsConstant;
	private boolean isLikelihood = true;
	
	private void init(byte gapPenalty) {
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
	
	public MatchScoring(byte gapPenalty) {
		init(gapPenalty);
	}
	
	public MatchScoring(String scoringJSON) {
		JsonParserFactory factory=JsonParserFactory.getInstance();
		JSONParser scrParser=factory.newJsonParser();
		Map scrMap = scrParser.parseJson(scoringJSON);
		
		String gapString = (String)scrMap.get(GAP_KEY);
		if (LINEAR_PENALTY_STR.equalsIgnoreCase(gapString)) {
			init(LINEAR_PENALTY);
		}
		else if (CONSTANT_PENALTY_STR.equalsIgnoreCase(gapString)) {
			init(CONSTANT_PENALTY);
		}
		else if (AFFINE_PENALTY_STR.equalsIgnoreCase(gapString)) {
			init(AFFINE_PENALTY);
		}
		else if (CONVEX_PENALTY_STR.equalsIgnoreCase(gapString)) {
			init(CONVEX_PENALTY);
		}
		else if (PROFILE_PENALTY_STR.equalsIgnoreCase(gapString)) {
			init(PROFILE_PENALTY);
		}
		else {
			throw new IllegalArgumentException("Unknown gap penalty: "+gapString);
		}

		if (scrMap.containsKey(GAPCONST_KEY)) {
			gapConst = Double.parseDouble((String)scrMap.get(GAPCONST_KEY));
		}
		if (scrMap.containsKey(GAPLIN_KEY)) {
			gapLinear = Double.parseDouble((String)scrMap.get(GAPLIN_KEY));
		}
		if (scrMap.containsKey(SUBCONST_KEY)) {
			subsConstant = Double.parseDouble((String)scrMap.get(SUBCONST_KEY));
		}
		if (scrMap.containsKey(SUBLIKLI_KEY)) {
			isLikelihood = "true".equalsIgnoreCase((String)scrMap.get(SUBLIKLI_KEY));
		}
		if (scrMap.containsKey(SUBMAP_KEY)) {
			setSubstitutionMatrix((Map)scrMap.get(SUBMAP_KEY));
		}
	}
	
	public String toJSON() {
		Map scrMap = new HashMap();
		switch(gapType) {
			case CONSTANT_PENALTY:
				scrMap.put(GAP_KEY, CONSTANT_PENALTY_STR);
				break;
			case LINEAR_PENALTY:
				scrMap.put(GAP_KEY, LINEAR_PENALTY_STR);
				break;
			case AFFINE_PENALTY:
				scrMap.put(GAP_KEY, AFFINE_PENALTY_STR);
				break;
			case CONVEX_PENALTY:
				scrMap.put(GAP_KEY, CONVEX_PENALTY_STR);
				break;
			case PROFILE_PENALTY:
				scrMap.put(GAP_KEY, PROFILE_PENALTY_STR);
				break;
		}
		scrMap.put(GAPCONST_KEY, ""+gapConst);
		scrMap.put(GAPLIN_KEY, ""+gapLinear);
		scrMap.put(SUBCONST_KEY, ""+subsConstant);
		scrMap.put(SUBLIKLI_KEY, isLikelihood ? "true" : "false");
		if (subMatrix != null) scrMap.put(SUBMAP_KEY, subMatrix);
		
		JsonGeneratorFactory genFactory=JsonGeneratorFactory.getInstance();
        JSONGenerator generator=genFactory.newJsonGenerator();
        String theAnswer = generator.generateJson(scrMap);
        return theAnswer.substring(1,theAnswer.length()-1);
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
		if ("X".equalsIgnoreCase(c1) || "X".equalsIgnoreCase(c2) && !subMatrix.containsKey("X")) {
			return Double.MAX_VALUE;
		}
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
