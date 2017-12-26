package circe.main;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.ReferenceGenome;

import com.javamex.classmexer.MemoryUtil;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableStream;

/**
 * Software to perform annotation independent detection of circular RNAs. 
 * @author Andrea Borbón and David Ayala Usma
 * REQUIRES A BAM FILE FROM A PAIRED-END RNA-Seq EXPERIMENT  
 */

public class CIRCE {
	
	//------------------------------------------------------------------------
	// Constants
	//------------------------------------------------------------------------
	
	public static final int CLIPPING_THRESHOLD = 19;
	public static final int MAX_DISTANCE_ALNS = 200000;
	public static final int MAX_ALLOWED_ALNS = 3;
	public static final int SPLICING_SIGNAL_TOLERANCE_WINDOW = 5;
	public static final int JUNCTION_BOUNDARY_COORDINATE_TOLERANCE_WINDOW = 20;
	
	//------------------------------------------------------------------------
	// Attributes
	//------------------------------------------------------------------------

	/**
	 * List of detected alignments.
	 */
	private HashMap<String, ArrayList<ReadAlignment>> alignments;
	
	/**
	 * Number of stored alignments so far
	 */
	private int storedAlignments;
	
	/**
	 * Reference genome of the organism
	 */
	private ReferenceGenome refGenome;
	
	/**
	 * 
	 */
	private ArrayList<CircRNA> predictedCircRNAs;
	
	//------------------------------------------------------------------------
	// Main methods
	//------------------------------------------------------------------------
	
	/**
	 * Main method to run the program
	 * @param args Array with one element, the path to the alignments file
	 * @throws Exception If the file can not be read
	 */
	public static void main(String[] args) throws Exception {
		CIRCE instance = new CIRCE();
		System.err.println("-------------------------------------------- CIRCE Output Log -------------------------------------------");
		System.err.println("[" + instance.getTimeStamp() + "]" + " Run started." );
		System.err.println("[" + instance.getTimeStamp() + "]" + " Loading reference genome." );
		instance.refGenome = new ReferenceGenome(args[1]);
		System.err.println("[" + instance.getTimeStamp() + "]" + " Scanning BAM file." );
		instance.storedAlignments = 0;
		instance.alignments = new HashMap<String, ArrayList<ReadAlignment>>();
		instance.processAlignmentsFile(args[0]);
		instance.printOutput();
		System.err.println("");
		System.err.println("[" + instance.getTimeStamp() + "]" + " Run finished." );
	}

	/**
	 * Read processor of the CIRCE algorithm
	 * @param filename Path to the BAM file to process
	 * @throws IOException If the file cannot be read
	 */
	public void processAlignmentsFile(String filename) throws IOException {		
		
		//Creating the alignment reader and writer files with HTSJDK.
		SamReader reader = null;
		SAMFileWriter BAMWriter = null;
		
		try 
		{
			//Creating the alignment file reader.
			reader = SamReaderFactory.makeDefault().open(new File(filename));
			
			//Creating the temporary output file and the alignment file writer.
			String tmpFilePath = System.getProperty("user.dir");
			File tmpBAMFile = new File(tmpFilePath + "tmp.bam");
			BAMWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, tmpBAMFile); 
			
			//Creating the iterator for the alignment file. 
			Iterator<SAMRecord> it = reader.iterator();
			
			//Counters for analyzed and saved alignments.
			int alignmentCounter = 1;
			int compliantAlignments = 0;
			
			//Reading the alignment file.			
			while (it.hasNext()) 
			{			
				//Recovers the next element.
				SAMRecord aln = it.next();
				
				//Looking for PCC signals in the alignments.
				if(!aln.getReadUnmappedFlag() && !aln.isSecondaryOrSupplementary() &&  aln.getCigarLength() > 1)
				{	
					int firstElementLength = aln.getCigar().getFirstCigarElement().getLength();
					int lastElementLength = aln.getCigar().getLastCigarElement().getLength();
					
					if((aln.getCigar().isLeftClipped() && firstElementLength >= CLIPPING_THRESHOLD) 
					   || (aln.getCigar().isRightClipped() && lastElementLength >= CLIPPING_THRESHOLD)) 
					{
						//If the alignment presents PCC signals, it is stored in the temporary file and the count increases.
						BAMWriter.addAlignment(aln);												
						compliantAlignments++;
					}
				}
				
				//Counts the number of processed alignments and reports the progress and memory usage per each 1.000.000 alignments.
				alignmentCounter++;
				
				if(alignmentCounter % 1000000 == 0)
				{
					DecimalFormat formatter = new DecimalFormat("###,###,###,###,###"); 
					String lineOutput = formatter.format(alignmentCounter);
					long memoryUsage = MemoryUtil.deepMemoryUsageOf(this) / 1048576;
					String storedAlns = formatter.format(compliantAlignments);
					System.err.println("[" + getTimeStamp() + "] " + lineOutput +" alignments processed - " + storedAlns + " alignments stored - " + memoryUsage + " MB of RAM used.");
				}
				
			}
			
			//Sets the final number of alignments that complied with the filters.
			storedAlignments = compliantAlignments;
		} 
		
		finally {
			
			if(reader != null || BAMWriter != null) 
			{
				reader.close();
				BAMWriter.close();
			}
		}
		
		//Removing the unique alignments that complied with the previous conditions.
		removeUniqueAlignments();
		
		//Sorting the ArrayLists of read alignments in the HashMap.
		sortingListsInAlignmentsMap();
		
		//Verifying distance and CIGAR criteria.
		filterAlignmentsByMaxDistanceAndCigar();
		
		//Filtering the junction read candidates by the location of their mate pairs.
		filterByMatePairLocation();
		
		//Filtering by splicing signals.
		filterBySplicingSignals();
		
		//Perform the circular RNA prediction with the filtered reads.
		predictCircularRNAs();
		
		//Final results
		recordNumberOfPredictedCircRNAs();
		
	}
	
	
	/**
	 * Method that remove those alignments that appear only once in the current alignments ArrayList.
	 */
	public void removeUniqueAlignments()
	{
		/**
		 * NUEVA ESTRATEGIA: UTILIZO UN DOBLE HASHING EN DONDE RECORRO EL BAM TEMPORAL PARA CONTAR EL NÚMERO DE ALINEAMIENTOS DE CADA LECTURA, 
		 * ALMACENANDO ESOS VALORES EN UNA TABLA DE HASH. DESPUÉS RECORRO EL ARCHIVO BAM TEMPORAL, RECUPERO EL NOMBRE DE LA LECTURA Y GUARDO
		 * EL ALINEAMIENTO EN LA TABLA DE ALINEAMIENTOS SI ESA LECTURA TIENE MÁS DE UN ALINEAMIENTO SEGÚN LA SEGUNDA TABLA. 
		 * ESTO LO HAGO CON LAS CLASES DE NGSEP.
		 */
		
		
		
		
		//Verbose response
		System.err.println("[" + getTimeStamp() + "]" + " BAM scanning finished. Starting the unique alignments filtering." );
		
		//Creation of the iterator through the whole HashMap and alignments tracking 
		Set<Map.Entry<String, ArrayList<ReadAlignment>>> entries = alignments.entrySet();
		Iterator<Map.Entry<String, ArrayList<ReadAlignment>>> iterator = entries.iterator();
		int counter = 0;
		
		//Removing the reads that only have one element
		while (iterator.hasNext()) 
		{
			Map.Entry<String, ArrayList<ReadAlignment>> currentRead = iterator.next();
			
			int numberOfReadAlignments = currentRead.getValue().size();
			counter += numberOfReadAlignments;
			
			if(numberOfReadAlignments == 1)
			{
				iterator.remove();
				counter -= 1;
			}
		}
		
		//Setting the value of stored alignments to the count made here
		storedAlignments = counter;
	}
	

	/**
	 * Performs the sorting of ArrayLists in the alignments HashMap.
	 */
	public void sortingListsInAlignmentsMap()
	{
		//Notification to user
		System.err.println("[" + getTimeStamp() + "]" + " Unique alignments filtering finished. Starting alignments sorting by coordinate." );
		
		//Recovering the reads stored in the collection
		Set<String> reads = alignments.keySet();
		
		for(String currentRead:reads) {
			
			ArrayList<ReadAlignment> currentReadAlignments = alignments.get(currentRead);
			sorterByCoordinate(currentReadAlignments);
		}
		
	}
	
	/**
	 * Filters the read in the HashMap by compliance with the maximum distance parameter and CIGAR structure.
	 */
	public void filterAlignmentsByMaxDistanceAndCigar()
	{
		//Notification to user
		System.err.println("[" + getTimeStamp() + "]" + " Alignments sorting by coordinate finished. Starting distance and CIGAR operators filtering." );
		
		//Creation of the iterator through the whole HashMap and alignments tracking 
		Set<Map.Entry<String, ArrayList<ReadAlignment>>> entries = alignments.entrySet();
		Iterator<Map.Entry<String, ArrayList<ReadAlignment>>> iterator = entries.iterator();
	
		while (iterator.hasNext()) 
		{
			//Recovering the alignments
			Map.Entry<String, ArrayList<ReadAlignment>> currentRead = iterator.next();
			ArrayList<ReadAlignment> currentReadAlignments = currentRead.getValue();
			
			//Distance and number of alignments calculation
			int numberAlignments = currentReadAlignments.size();
			ReadAlignment firstAln = currentReadAlignments.get(0);
			ReadAlignment lastAln = currentReadAlignments.get(numberAlignments - 1);
			int distanceFirstLastAlns = 0;
						
			//Strand sense verification
			if(firstAln.isPositiveStrand() && lastAln.isPositiveStrand())
			{
				distanceFirstLastAlns = lastAln.getFirst() - firstAln.getFirst();				
			}
			
			else if(firstAln.isNegativeStrand() && lastAln.isNegativeStrand())
			{
				distanceFirstLastAlns = lastAln.getLast() - firstAln.getLast();
			}

			//Distance and alignment number filtering
			if(numberAlignments > MAX_ALLOWED_ALNS || distanceFirstLastAlns > MAX_DISTANCE_ALNS)
			{
				iterator.remove();
			}
			
			//CIGAR verification
			else {
				int firstAlnNumOperators = firstAln.getNumCigarItems();
				int lastAlnNumOperators = lastAln.getNumCigarItems();
				
				int firstAlnCigarOperator = firstAln.getCigarItemOperator(0);
				int lastAlnCigarOperator = lastAln.getCigarItemOperator(lastAln.getNumCigarItems() - 1);
				
				//If the read does not comply with the conditions that the Leftmost alignment CIGAR == H/S AND Rightmost alignment CIGAR == H/S, it must be removed. 
				if((firstAlnNumOperators > 2 || lastAlnNumOperators > 2) || 
					!((firstAlnCigarOperator == ReadAlignment.ALIGNMENT_HARDCLIP || firstAlnCigarOperator == ReadAlignment.ALIGNMENT_SKIPFROMREAD) 
					&& (lastAlnCigarOperator == ReadAlignment.ALIGNMENT_HARDCLIP || lastAlnCigarOperator == ReadAlignment.ALIGNMENT_SKIPFROMREAD)))
				{
					iterator.remove();
				}
			}
		}
	}
	
	/**
	 * Filters the reads in the HashMap 
	 */
	public void filterByMatePairLocation()
	{
		//Notification to user
		System.err.println("[" + getTimeStamp() + "]" + " Distance and CIGAR operators filtering finished. Starting the verification of read mates in the experiment." );
		
		//Creation of the iterator through the whole HashMap 
		Set<Map.Entry<String, ArrayList<ReadAlignment>>> entries = alignments.entrySet();
		Iterator<Map.Entry<String, ArrayList<ReadAlignment>>> iterator = entries.iterator();
	
		while (iterator.hasNext()) 
		{
			//Recovering the alignments of the current read.
			Map.Entry<String, ArrayList<ReadAlignment>> currentRead = iterator.next();
			ArrayList<ReadAlignment> currentReadAlignments = currentRead.getValue();
			int numberOfAlignments = currentReadAlignments.size();
			ReadAlignment firstAlignment = currentReadAlignments.get(0);
			ReadAlignment lastAlignment = currentReadAlignments.get(numberOfAlignments - 1);
			
			//Verification of the mate position
			if(!(firstAlignment.isPaired() && lastAlignment.isPaired() && firstAlignment.getFirst() < firstAlignment.getMateFirst() && lastAlignment.getFirst() > lastAlignment.getMateFirst()))
			{
				iterator.remove();
			}
		}
	}
	
	public void filterBySplicingSignals()
	{
		//Notification to user
		System.err.println("[" + getTimeStamp() + "]" + " Verification of read mates in the experiment finished. Splicing signal filtering started." );
		
		//Creation of the iterator through the whole HashMap 
		Set<Map.Entry<String, ArrayList<ReadAlignment>>> entries = alignments.entrySet();
		Iterator<Map.Entry<String, ArrayList<ReadAlignment>>> iterator = entries.iterator();
	
		while (iterator.hasNext()) 
		{
			//Recovering the alignments of the current read.
			Map.Entry<String, ArrayList<ReadAlignment>> currentRead = iterator.next();
			ArrayList<ReadAlignment> currentReadAlignments = currentRead.getValue();
			int numberOfAlignments = currentReadAlignments.size();
			ReadAlignment firstAln = currentReadAlignments.get(0);
			ReadAlignment lastAln = currentReadAlignments.get(numberOfAlignments - 1);
			
			int firstAlnLastPos = firstAln.getLast();
			int lastAlnFirstPos = lastAln.getFirst();
			
			if(firstAln.getSequenceName().equals(lastAln.getSequenceName()))
			{
				String sequenceName = firstAln.getSequenceName();
				
				StringBuilder acceptorSite = new StringBuilder(); 
				StringBuilder donorSite = new StringBuilder();
				
				//Recovering the sequence for the splicing tolerance window
				for (int i = 1; i <= SPLICING_SIGNAL_TOLERANCE_WINDOW; i++) 
				{
					donorSite.append(refGenome.getReferenceBase(sequenceName, firstAlnLastPos + i));
					acceptorSite.append(refGenome.getReferenceBase(sequenceName, (lastAlnFirstPos - SPLICING_SIGNAL_TOLERANCE_WINDOW - 1) + i));
				}
					
				String acceptorSiteWindow = acceptorSite.toString();
				String donorSiteWindow = donorSite.toString();
				
				//Verify strand-specific splicing signals in the reads
				if(firstAln.isPositiveStrand() && lastAln.isPositiveStrand())
				{
					if(!(acceptorSiteWindow.contains("AG") && donorSiteWindow.contains("GT")))
					{
						iterator.remove();	
					}
				}
				
				else if (firstAln.isNegativeStrand() && lastAln.isNegativeStrand())
				{
					if(!(acceptorSiteWindow.contains("AC") && donorSiteWindow.contains("CT")))
					{
						iterator.remove();
					}
				}
				
				else
				{
					iterator.remove();
				}
			}

			else
			{
				iterator.remove();
			}
		}	
	}
	
	public void predictCircularRNAs()
	{
		//Notification to user
		System.err.println("[" + getTimeStamp() + "]" + " Splicing signal filtering finished. circRNA prediction started." );
		
		//Creating the array of CircRNAs and preparing the control structure of the visited reads
		predictedCircRNAs = new ArrayList<CircRNA>();
		Set<String> countedReads = new HashSet<String>();
		ArrayList<String> allAlignments = new ArrayList<String>(alignments.keySet());
		int circIDs = 0;
		
		for(int i = 0; i < allAlignments.size(); i++)
		{	
			//Skipping this value if this read has been already counted.
			String currentRead = allAlignments.get(i);
			if(countedReads.contains(currentRead))
			{
				continue;
			}
			
			//Features of the new circRNA
			String currentReferenceSequence = alignments.get(currentRead).get(0).getSequenceName();
			int currentReadAlignments = alignments.get(currentRead).size();
			int currentStartCoordinate = alignments.get(currentRead).get(0).getFirst();
			int currentEndCoordinate = alignments.get(currentRead).get(currentReadAlignments - 1).getLast();
			char currentCodingStrand = '+';
			if(alignments.get(currentRead).get(0).isNegativeStrand())
			{
				currentCodingStrand = '-';
			}
			
			//Creation of the circRNA entry
			CircRNA newCircRNA = new CircRNA("circRNA" + circIDs, currentReferenceSequence, currentStartCoordinate, currentEndCoordinate, currentCodingStrand, 1, currentRead);

			
			//Looking for the supporting reads and removing them from the 
			for(int j = i + 1; j < allAlignments.size(); j++)
			{
				//Once again skipping this value if this read has been already counted.
				String inspectedRead = allAlignments.get(j);
				if(countedReads.contains(inspectedRead))
				{
					continue;
				}
				
				//Features of the inspected read
				String inspectedReferenceSequence = alignments.get(inspectedRead).get(0).getSequenceName();
				int inspectedReadAlignments = alignments.get(inspectedRead).size();
				int inspectedStartCoordinate = alignments.get(inspectedRead).get(0).getFirst();
				int inspectedEndCoordinate = alignments.get(inspectedRead).get(inspectedReadAlignments - 1).getLast();
				char inspectedCodingStrand = '+';
				if(alignments.get(inspectedRead).get(0).isNegativeStrand())
				{
					inspectedCodingStrand = '-';
				}
				
				//Verification if this read supports the newly created circRNA: Within the same boundary window coordinates, in the same sequence, and in the same strand. 
				boolean supportingRead = Math.abs(currentStartCoordinate - inspectedStartCoordinate) <= JUNCTION_BOUNDARY_COORDINATE_TOLERANCE_WINDOW && Math.abs(currentEndCoordinate - inspectedEndCoordinate) <= JUNCTION_BOUNDARY_COORDINATE_TOLERANCE_WINDOW;				
			
				if(!(supportingRead && inspectedReferenceSequence.equals(currentReferenceSequence) && inspectedCodingStrand == currentCodingStrand))
				{
					continue;
				}
				
				else
				{
					//Marking the read as visited
					countedReads.add(inspectedRead);
					
					//Replacing the circRNA coordinates if the inspected leftmost is less than the current leftmost and the inspected rightmost is greater than the current rightmost
					if(inspectedStartCoordinate < currentStartCoordinate)
					{
						newCircRNA.setStartCoordinate(inspectedStartCoordinate);
					}
					
					if(inspectedEndCoordinate > currentEndCoordinate)
					{
						newCircRNA.setEndCoordinate(inspectedEndCoordinate);
					}
					
					//Adding a support read to the count
					newCircRNA.setNumberJunctionReadsSupport(newCircRNA.getNumberJunctionReadsSupport() + 1);
					
					//Adding a support read name to the count
					newCircRNA.setNameOfSupportingJunctionReads(inspectedRead);
					
				}
				
			}
			
			//Adding the circRNA to the list.
			predictedCircRNAs.add(newCircRNA);
			circIDs++;
		}
		
	}
	
	/**
	 * Refining the splicing junction coordinates using splicing signals 
	 */
	public void refiningCoordinates() 
	{
		//Notification to user
		System.err.println("[" + getTimeStamp() + "]" + " circRNA prediction finished. Coordinate refining started." );
		
		for (CircRNA circRNA : predictedCircRNAs) 
		{
			int leftCoordinate = circRNA.getStartCoordinate();
			int rightCoordinate = circRNA.getEndCoordinate();
			char codingStrand = circRNA.getCodingStrand();
			String sequenceName = circRNA.getNameReferenceSequence();
			
			StringBuilder leftSite = new StringBuilder();
			StringBuilder rightSite = new StringBuilder();
			
			//Recovering the sequence for the splicing tolerance window
			for (int i = 1; i <= SPLICING_SIGNAL_TOLERANCE_WINDOW; i++) 
			{
				rightSite.append(refGenome.getReferenceBase(sequenceName, rightCoordinate + i));
				leftSite.append(refGenome.getReferenceBase(sequenceName, (leftCoordinate - SPLICING_SIGNAL_TOLERANCE_WINDOW - 1) + i));
			}
				
			String leftSiteWindow = leftSite.toString();
			String rightSiteWindow = rightSite.toString();
			
			//Verify strand-specific splicing signals in the reads
			if(codingStrand == '+')
			{
				int coordinateSignalLeft = leftCoordinate - leftSiteWindow.indexOf("AG") + 2;
				int coordinateSignalRight = rightCoordinate + rightSiteWindow.indexOf("GT") - 2;
				circRNA.setStartCoordinate(coordinateSignalLeft);
				circRNA.setEndCoordinate(coordinateSignalRight);
			}
			
			else
			{
				int coordinateSignalLeft = leftCoordinate - leftSiteWindow.indexOf("AC") + 2;
				int coordinateSignalRight = rightCoordinate + rightSiteWindow.indexOf("CT") - 2;
				circRNA.setStartCoordinate(coordinateSignalLeft);
				circRNA.setEndCoordinate(coordinateSignalRight);
			}
		}
	}
	
	//------------------------------------------------------------------------
	// Helper methods
	//------------------------------------------------------------------------
	
	/**
	 * Prints the number of predicted circRNAs to standard output in the console
	 */
	public void recordNumberOfPredictedCircRNAs()
	{
		System.out.println("circRNA ID"+ "\t" + "Contig/chromosome" + "\t" + "Start coordinate" + "\t" + "End coordinate" + "\t" + "Coding strand" + "\t" + "Number of supporting junction reads" + "\t" + "Name of the supporting junction reads");
		for (CircRNA circRNA : predictedCircRNAs) 
		{
			System.out.println(circRNA.getCircRNAIdentifier() + "\t" + circRNA.getNameReferenceSequence() + "\t" + circRNA.getStartCoordinate() + "\t" + circRNA.getEndCoordinate() + "\t" + circRNA.getCodingStrand() + "\t" + circRNA.getNumberJunctionReadsSupport() + "\t" + circRNA.getNameOfSupportingJunctionReads());
		}
	}
	
	
	/**
	 * Sorts an ArrayList of ReadAlignments by coordinate of the first base in the alignment.
	 * @param alignmentsList - List to be sorted.
	 */
	public void sorterByCoordinate(ArrayList<ReadAlignment> alignmentsList)
	{
		//Creating a new comparator with a lambda function
		Comparator<ReadAlignment> compareByCoordinate = (ReadAlignment read1, ReadAlignment read2) -> Integer.compare(read1.getFirst(), read2.getFirst());
		
		//Sorting the array
		alignmentsList.sort(compareByCoordinate);
	}
	
		
	/**
	 * Prints the timestamp.
	 * @return String - Timestamp.
	 */
	public String getTimeStamp ()
	{
		String formattedDate = null;
		Date date = new Date();
		SimpleDateFormat sdf = new SimpleDateFormat("dd/MM/yyyy h:mm:ss a");
		formattedDate = sdf.format(date);
		return formattedDate;
	}
	
	/**
	 * Returns the complete memory size of the query object in megabytes.
	 * @param object Object to be measured
	 * @return long - Size of the measured object.
	 */
	public long getMemoryUsageOf(Object object)
	{
		long memoryUsage = MemoryUtil.deepMemoryUsageOf(object) / 1048576;
		return memoryUsage;
	}
	
	
	/**
	 * Prints the final statistics to standard outsput.
	 */
	public void printOutput() 
	{
		DecimalFormat formatter = new DecimalFormat("###,###,###,###,###"); 
		String predictedRNAs = formatter.format(predictedCircRNAs.size());
		System.err.print("[" + getTimeStamp() + "] " +  predictedRNAs + " predicted circRNAs.");

	}

}

