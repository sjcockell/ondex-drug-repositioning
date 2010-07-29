package net.sourceforge.ondex.transformer.version1;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.regex.Pattern;

import java.io.BufferedReader;
import java.io.FileReader;

import net.sourceforge.ondex.core.AttributeName;
import net.sourceforge.ondex.core.CV;
import net.sourceforge.ondex.core.ConceptClass;
import net.sourceforge.ondex.core.EvidenceType;
import net.sourceforge.ondex.core.ONDEXConcept;
import net.sourceforge.ondex.core.ONDEXGraph;
import net.sourceforge.ondex.core.ONDEXRelation;
import net.sourceforge.ondex.core.RelationType;
import net.sourceforge.ondex.core.ConceptAccession;
import net.sourceforge.ondex.event.type.AttributeNameMissingEvent;
import net.sourceforge.ondex.event.type.CVMissingEvent;
import net.sourceforge.ondex.event.type.ConceptClassMissingEvent;
import net.sourceforge.ondex.event.type.EvidenceTypeMissingEvent;
import net.sourceforge.ondex.event.type.RelationTypeMissingEvent;
import net.sourceforge.ondex.transformer.ONDEXTransformer;
import net.sourceforge.ondex.args.ArgumentDefinition;
import net.sourceforge.ondex.args.StringArgumentDefinition;

/**
 * Makes transformations necessary for converting dataset v1 into dataset v2
 * @author sjcockell
 * 
 */

public class Transformer extends ONDEXTransformer implements ArgumentNames {
    
	/**
	 * Returns name of this transformer.
	 *
	 * @return String
	 */
	public String getName() {
		return "ETherapeutics Version Updater";
	}
    public String getId() {
        return "version1";
    }
    /**
	 * Returns version of this transformer.
	 *
	 * @return String
	 */
	public String getVersion() {
	 	return "0.1.3";
	}

	/**
	 * Returns arguments required by this transformer.
	 *
	 * @return ArgumentDefinition<?>[]
	 */
	@Override
    public ArgumentDefinition<?>[] getArgumentDefinitions() {
        /*
        new ArgumentDefinition<?>[] {
        new StringArgumentDefinition(params)
        }
         */
        return new ArgumentDefinition<?>[]{
			new StringArgumentDefinition(INTERACTION_FILE_ARG, INTERACTION_FILE_ARG_DESC, true, null, false),
			new StringArgumentDefinition(MAPPING_FILE_ARG, MAPPING_FILE_ARG_DESC, true, null, false),
		};
    }


    /**
     * Does not require index ondex graph.
     *
     * @return false
     */
    public boolean requiresIndexedGraph() {
		return false;
	}

    public String[] requiresValidators() {
	    return new String[0];
	}
	
	public void start() throws Exception {
		String interactionFileName = convertToString(INTERACTION_FILE_ARG);
		String mappingFileName = convertToString(MAPPING_FILE_ARG);
		Map<String, Float> scores = getInteractionScores(interactionFileName);
		Map<String, String> mapping = getHPRDMapping(mappingFileName);
		CV eTheraCV = requireCV("eTherapeutics");
		AttributeName scoreGDS = requireAttributeName("G-Sesame Score");
		//all the Protein concepts
		Set<ONDEXConcept> conceptview = graph.getConceptsOfConceptClass(graph.getMetaData().getConceptClass("Protein"));
		//all the interaction relations
		Set<ONDEXRelation> relationview = graph.getRelationsOfRelationType(graph.getMetaData().getRelationType("it_wi"));
		for (ONDEXConcept c : conceptview) {
			//check the HPRD ID on each, and add eThera ID where appropriate (need a Hash of HPRD vs eThera IDs)
			String accession = getHPRDAccession(c);
			if (accession != null) {
				String eTheraID = mapping.get(accession);
				if (eTheraID != null) {
					c.createConceptAccession(eTheraID, eTheraCV, false);
				}
			}
		}
		for (ONDEXRelation r : relationview) {
			//need to add G-Sesame scores where appropriate - how?
			ONDEXConcept c1 = r.getFromConcept();
			ONDEXConcept c2 = r.getToConcept();
			String acc1 = getHPRDAccession(c1);
			String acc2 = getHPRDAccession(c2);
			//use 2 accessions to get G-Sesame score
			//getScore(acc1, acc2);
			Float s = scores.get(acc1+":"+acc2);
			//r.createRelationGDS();
            if (s != null) {
				r.createGDS(scoreGDS, s, false);
			}
		}
	}
	
	public String getHPRDAccession(ONDEXConcept c) {
		String acc = null;
		Set<ConceptAccession> ci = c.getConceptAccessions();
		for (ConceptAccession a : ci) {
			CV conceptCV = a.getElementOf();
			if (conceptCV.getId() == "hprd_id") {
				acc = a.getAccession();
			}
		}
		return acc;
	}

    public String convertToString(String arg) throws Exception { 
        Object o = ta.getUniqueValue(arg);
        if (o == null) {
            throw new IllegalArgumentException("Argument "+arg+" is required, but was not supplied");
        }
        String value = o.toString();
        return value;
    }


	public Map<String, Float> getInteractionScores(String intFile) throws Exception {
		Map<String, Float> scoreMap = new Hashtable<String, Float>();
		//open file with interaction-score info
		BufferedReader br = new BufferedReader(new FileReader(intFile));
		//read through file, creating int:int-score pairs
		String line = null;
        while ((line = br.readLine()) != null) {
			String[] columns = line.trim().split("\t");
			String hprdID1 = columns[0];
			String hprdID2 = columns[1];
			String score = columns[2];
			String ids = hprdID1+":"+hprdID2;
			Float fscore = null;
			if (!score.equals("None")) {
				fscore = Float.valueOf(score);
				scoreMap.put(ids, fscore);
			}
		}
		//return resulting Map
		return scoreMap;
	}

	public Map<String, String> getHPRDMapping(String mapFile) throws Exception {
		Map<String, String> idMap = new Hashtable<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(mapFile));
		String line = br.readLine(); //gets rid of first line - header.
		while ((line = br.readLine()) != null) {
			String[] columns = line.trim().split("\t");
			if (line.startsWith("hw_etdb")) {
				//there are some blank lines in the file, this means only info lines should be processed
				String hprdID = columns[1];
				String eTheraID = columns[0];
				if (hprdID != null && eTheraID != null) {
					idMap.put(hprdID, eTheraID);
				}
			}
		}
		return idMap;
	}
}

//code at the moment is ALL from a uniprot transformer...
/*public class Transformer {
	private ONDEXGraph graph;
	private ParserArguments pa;

	private ConceptClass ccProtein;

	private ConceptClass ccPublication;

	private ConceptClass ccEC;
	
	private ConceptClass ccDisease;

	private AttributeName attTaxId;

	private AttributeName attSequence;

	private AttributeName attTitle;
	
	private AttributeName attYear;
	
	private AttributeName attJournal;
	
	private AttributeName pubType;
	
	private HashMap<String, CV> cvs = new HashMap<String, CV>();

	private CV cvUniProt;
	private CV cvPubMed;
	private CV cvEC;
	private CV cvOMIM; 
	private CV cvDOI; 

	private EvidenceType etAutomaticallyCurated;
	private EvidenceType etManuallyCurated;
	private EvidenceType ev;
	
	private RelationType rtPublishedIn;
	private RelationType rtCatC;
	private RelationType rtInvIn;
	
	private boolean addContextInformation = false;
		
	private HashMap<String, Integer> ecToConcept = new HashMap<String, Integer>();
	private HashMap<String, Integer> omimToConcept = new HashMap<String, Integer>();
	private HashMap<String, Integer> pubmedToConcept = new HashMap<String, Integer>();
	
	public static HashSet<String> unknownCVs = new HashSet<String>();

	
	private Set<String> umambigProtinMetaData;
	
	
	public Transformer(ONDEXGraph graph, ParserArguments pa, boolean addContextInformation) {
		
	}
	
	
	public void transform() {

	}
*/	
