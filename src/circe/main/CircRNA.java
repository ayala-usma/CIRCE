package circe.main;

public class CircRNA 
{
	//------------------------------------------------------------------------
	// Attributes
	//------------------------------------------------------------------------

	/**
	 * Identifiers of the circular RNA identified
	 */
	private String circRNAIdentifier;
	
	/**
	 * The name of the reference sequence of the circRNA
	 */
	private String nameReferenceSequence;
	
	/**
	 * Leftmost coordinate of the junction according to the reference genome
	 */
	private int startCoordinate;
	
	/**
	 * Rightmost coordinate of the junction according to the reference genome
	 */
	private int endCoordinate;
	
	/**
	 * Strand that codes the circRNA
	 */
	private char codingStrand;
	
	/**
	 * Number of junction reads that support the circRNA.
	 */
	private int numberJunctionReadsSupport;
	
	/**
	 * Name of the junction reads that support this circRNA
	 */
	private String nameOfSupportingJunctionReads;
	
	
	//------------------------------------------------------------------------
	// Constructor
	//------------------------------------------------------------------------
	
	public CircRNA(String circRNAIdentifier, String nameReferenceSequence, int startCoordinate, int endCoordinate,
			char codingStrand, int numberJunctionReadsSupport, String nameOfSupportingJunctionReads) 
	{
		this.circRNAIdentifier = circRNAIdentifier;
		this.nameReferenceSequence = nameReferenceSequence;
		this.startCoordinate = startCoordinate;
		this.endCoordinate = endCoordinate;
		this.codingStrand = codingStrand;
		this.numberJunctionReadsSupport = numberJunctionReadsSupport;
		this.nameOfSupportingJunctionReads = nameOfSupportingJunctionReads;
	}


	//------------------------------------------------------------------------
	// Methods
	//------------------------------------------------------------------------
	
	/**
	 * Returns the circRNA identifier
	 * @return String - circRNA identifier
	 */
	public String getCircRNAIdentifier() 
	{
		return circRNAIdentifier;
	}


	/**
	 * Sets the circRNA identifier
	 * @param circRNAIdentifier
	 */
	public void setCircRNAIdentifier(String circRNAIdentifier) 
	{
		this.circRNAIdentifier = circRNAIdentifier;
	}

	/**
	 * Returns the reference sequence of this circRNA
	 * @return String - Reference sequence of this circRNA
	 */
	public String getNameReferenceSequence() 
	{
		return nameReferenceSequence;
	}

	/**
	 * Sets the reference sequence of this circRNA
	 * @param nameReferenceSequence
	 */
	public void setNameReferenceSequence(String nameReferenceSequence) 
	{
		this.nameReferenceSequence = nameReferenceSequence;
	}

	/**
	 * Returns the start coordinate of this circRNA
	 * @return int - Start coordinate of this circRNA
	 */
	public int getStartCoordinate() 
	{
		return startCoordinate;
	}

	/**
	 * Sets the start coordinate of this circRNA
	 * @param startCoordinate
	 */
	public void setStartCoordinate(int startCoordinate) 
	{
		this.startCoordinate = startCoordinate;
	}

	/**
	 * Returns the end coordinate of this circRNA
	 * @return int - End coordinate of this circRNA
	 */
	public int getEndCoordinate() 
	{
		return endCoordinate;
	}

	/**
	 * Sets the end coordinate of this circRNA
	 * @param endCoordinate
	 */
	public void setEndCoordinate(int endCoordinate) 
	{
		this.endCoordinate = endCoordinate;
	}

	/**
	 * Returns the coding strand of this circRNA
	 * @return - char Coding strand of this circRNA
	 */
	public char getCodingStrand() 
	{
		return codingStrand;
	}

	/**
	 * Sets the coding strand of this circRNA
	 * @param codingStrand
	 */
	public void setCodingStrand(char codingStrand) 
	{
		this.codingStrand = codingStrand;
	}

	/**
	 * Returns the number of junction reads that support this circRNA
	 * @return int - Number of junction reads that support this circRNA.
	 */
	public int getNumberJunctionReadsSupport() 
	{
		return numberJunctionReadsSupport;
	}

	/**
	 * Sets the number of junction reads that support this circRNA
	 * @param numberJunctionReadsSupport
	 */
	public void setNumberJunctionReadsSupport(int numberJunctionReadsSupport) 
	{
		this.numberJunctionReadsSupport = numberJunctionReadsSupport;
	}
	
	/**
	 * Returns the name of the reads that support this circRNA
	 * @return String - Name of the reads that support this circRNA
	 */
	public String getNameOfSupportingJunctionReads()
	{
		return nameOfSupportingJunctionReads;
	}
	
	public void setNameOfSupportingJunctionReads(String newReadName)
	{
		this.nameOfSupportingJunctionReads = this.nameOfSupportingJunctionReads + "," + newReadName;
	}
	
}
