package dea;

public enum RankingType {
	STANDARD ("Thus if A ranks ahead of B and C (which compare equal) which are both ranked ahead of D," +
			"then A gets ranking number 1 (first), B gets ranking number 2 (joint second)," +
			"C also gets ranking number 2 (joint second)and D gets ranking number 4 (fourth)." +
			"In this case, nobody would get ranking number 3 (third) and that would be left as a gap."),
	DENSE ("Thus if A ranks ahead of B and C (which compare equal) which are both ranked ahead of D," +
			"then A gets ranking number 1 (first), B gets ranking number 2 (joint second)," +
			"C also gets ranking number 2 (joint second) and D gets ranking number 3 (third)");
	
	
	private String rankingDescription;
	RankingType(String Desc) {
		this.rankingDescription = Desc;
	}
	
	//Method to access _ModelType
	public String GetStatusDescription()
	{
		return rankingDescription;
	}
	
}
