package net.sourceforge.ondex.parser.symatlas;

import net.sourceforge.ondex.parser.ONDEXParser;
import net.sourceforge.ondex.args.ArgumentDefinition;
import net.sourceforge.ondex.args.StringArgumentDefinition;
import net.sourceforge.ondex.args.FileArgumentDefinition;
import net.sourceforge.ondex.core.ONDEXConcept;
import net.sourceforge.ondex.core.ONDEXRelation;
import net.sourceforge.ondex.core.ConceptAccession;
import net.sourceforge.ondex.core.ConceptClass;
import net.sourceforge.ondex.core.RelationType;
import net.sourceforge.ondex.core.CV;
import net.sourceforge.ondex.core.EvidenceType;
import net.sourceforge.ondex.core.AttributeName;

import net.sourceforge.ondex.exception.type.MetaDataMissingException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Hashtable;
import java.util.Map;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;

public class Parser extends ONDEXParser implements ArgumentNames {
    /* 
    Author - sjcockell
    Date - 20090630
    Version - 0.0.1
     */

    private ConceptClass proteinCC;
    private ConceptClass affyCC;
    private RelationType proteinAffyRelation;
    private CV genericCV;
    private EvidenceType genericET;
    private AttributeName scoreGDS;
    private AttributeName expressionGDS;
	private CV uniprotCV;
    private CV hprdCV;
    private CV hugoCV;
    private CV affyCV;
    private ONDEXConcept MatchedConcept = null;
    private Map<String, ONDEXConcept> proteins = new Hashtable<String, ONDEXConcept>();
    private Map<String, String> uniprotInfo = new Hashtable<String, String>();
    private Map<String, String> hprdInfo = new Hashtable<String, String>();
    private Map<String, String> hugoInfo = new Hashtable<String, String>();
    private Map<String, String> gsInfo = new Hashtable<String, String>();
    private Map<String, String> gnInfo = new Hashtable<String, String>();
    private Map<String, String> affyMapping = new Hashtable<String, String>(); 
	public boolean readsDirectory() {
        return true;
    }

    public boolean readsFile() {
        return false;
    }

    public ArgumentDefinition<?>[] getArgumentDefinitions() {
        /*
        new ArgumentDefinition<?>[] {
        new StringArgumentDefinition(params)
        }
         */
		return new ArgumentDefinition<?>[]{
			new StringArgumentDefinition(TISSUE_TYPE_ARG, TISSUE_TYPE_ARG_DESC, true, null, false),
            new FileArgumentDefinition(FileArgumentDefinition.INPUT_DIR, "Directory with input files", true, true, true, false)
		};
    }

    public String getName() {
        return "SymAtlas";
    }
    public String getId() {
        return "symatlas";
    }

    public String getVersion() {
        return "0.1.3";
    }

    public String[] requiresValidators() {
        return new String[0];
    }

    public boolean requiresIndexedGraph() {
        return false;
    }

    public void start() throws Exception {
        fetchMetadata();
		int concept_id = getLastIdForConcept();
		String[] tissue_type_list = parseArgument(TISSUE_TYPE_ARG);
        //get directory name
        String dirname = (String)pa.getUniqueValue(FileArgumentDefinition.INPUT_DIR);
		//affy annotation needs to be affyid:uniprotid
        affyMapping = readAffyAnnotation(dirname + "/uniprot_affy");
		Set<ONDEXConcept> conceptview = graph.getConceptsOfConceptClass(graph.getMetaData().getConceptClass("Protein"));
		Map<String, Integer> conceptMapping = setConceptMapping(conceptview);
    	Map<String, ONDEXConcept> probes = new Hashtable<String, ONDEXConcept>();
        for (String tissue_type:tissue_type_list) {
			//get the expression data for this tissue
			Map<String, Float> tissueExpressionData = readExpressionData(dirname + "/expression", tissue_type);
			//create Attribute for tissue expression data (need to make sure they're all in metadata)
			String attName = tissue_type + " EXPRESSION VALUE";
			AttributeName tissueAttribute = requireAttributeName(attName);
			//for every expression value, set up the concepts (or retrieve them
			//from the graph if they already exist)
			Iterator it = tissueExpressionData.keySet().iterator();
			while (it.hasNext()) {
				String affyid = (String)it.next();
				Float exp = tissueExpressionData.get(affyid);
				String uniprot_mapped = affyMapping.get(affyid);
                String[] mapped_array = null;
                if (uniprot_mapped != null) {
                    mapped_array = uniprot_mapped.split(":");
                }
                else {
                    mapped_array = new String[0];
                }
                for (String mapped:mapped_array) {
				    Integer existing_concept_id = null;
                    try {
                        existing_concept_id = conceptMapping.get(mapped);
			    	}
                    catch (NullPointerException npe) {
                        existing_concept_id = null;
                    }
                    if (existing_concept_id != null) {
					    //initially just link probe directly to protein - for quickest possible test...
					    ONDEXConcept c1 = graph.getConcept(existing_concept_id);
					    //need to check against a list of probes before doing this... (may be 2nd+ tissue type)
					    ONDEXConcept c2 = probes.get(affyid);
					    if (c2 == null) {
						    c2 = graph.getFactory().createConcept(++concept_id+"", affyCV, affyCC, genericET);
						    c2.createConceptName(affyid, true);
						    c2.createConceptAccession(affyid, affyCV, false);
						    c2.createGDS(tissueAttribute, exp, false);
						    probes.put(affyid, c2);
				    	}
					    else {
						    c2.createGDS(tissueAttribute, exp, false);
					    }
                        ONDEXRelation r1 = graph.getFactory().createRelation(c1, c2, proteinAffyRelation, genericET);
				    }
                }
			}
		}
    }

	private String[] parseArgument(String arg) throws Exception {
		Object o = pa.getUniqueValue(arg);
		if (o == null) {
			throw new IllegalArgumentException("Argument "+arg+" is required, but was not supplied");
		}
		String value = o.toString();
		ArrayList tissue_list = new ArrayList();
		if (value.equalsIgnoreCase("all")) {
			tissue_list = getDefaultTissueList();
		}
		else {
			String[] tokens = value.split(",");
			for (String token:tokens) {
				//strip spaces from start/end
				token = token.trim();
				tissue_list.add(token);
			}
		}
		String[] list = new String[tissue_list.size()];
        int i = 0;
        for (Object ob : tissue_list) { 
            list[i] = (String) ob;
            i++;
        }
		return list;
	}

	private ArrayList getDefaultTissueList() {
		String[] list = {
			"ColorectalAdenocarcinoma","WHOLEBLOOD","BM-CD33+Myeloid",
			"PB-CD14+Monocytes","PB-BDCA4+Dentritic_Cells","PB-CD56+NKCells",
			"PB-CD4+Tcells","PB-CD8+Tcells","PB-CD19+Bcells","BM-CD105+Endothelial",
			"BM-CD34+","leukemialymphoblastic(molt4)","721_B_lymphoblasts",
			"lymphomaburkittsRaji","leukemiapromyelocytic(hl60)","lymphomaburkittsDaudi",
			"leukemiachronicmyelogenous(k562)","thymus","Tonsil","lymphnode",
			"fetalliver","BM-CD71+EarlyErythroid","bonemarrow","TemporalLobe",
			"globuspallidus","CerebellumPeduncles","cerebellum","caudatenucleus",
			"WholeBrain","ParietalLobe","MedullaOblongata","Amygdala",
			"PrefrontalCortex","OccipitalLobe","Hypothalamus","Thalamus",
			"subthalamicnucleus","CingulateCortex","Pons","spinalcord",
			"fetalbrain","adrenalgland","Lung","Heart","Liver","kidney","Prostate",
			"Uterus","Thyroid","fetalThyroid","fetallung","PLACENTA",
			"CardiacMyocytes","SmoothMuscle","bronchialepithelialcells","ADIPOCYTE",
			"Pancreas","PancreaticIslets","testis","TestisLeydigCell",
			"TestisGermCell","TestisInterstitial","TestisSeminiferousTubule",
			"salivarygland","trachea","AdrenalCortex","Ovary","Appendix","skin",
			"ciliaryganglion","TrigeminalGanglion","atrioventricularnode","DRG",
			"SuperiorCervicalGanglion","SkeletalMuscle","UterusCorpus","TONGUE",
			"OlfactoryBulb","Pituitary"
		};
		return new ArrayList(Arrays.asList(list));
	}

    private void fetchMetadata() throws MetaDataMissingException {
        proteinCC = requireConceptClass("Protein");
        affyCC = requireConceptClass("AffymetrixProbe");
		proteinAffyRelation = requireRelationType("binds_to_encoding_mrna");
        genericCV = requireCV("eTherapeutics");
        uniprotCV = requireCV("UNIPROTKB");
        hprdCV = requireCV("HPRD");
        hugoCV = requireCV("HGNC");
		affyCV = requireCV("AFFYMETRIX");
        genericET = requireEvidenceType("IMPD");
        scoreGDS = requireAttributeName("G-Sesame Score");
		expressionGDS = requireAttributeName("Average Expression");
    }

	private Map<String, Integer> setConceptMapping(Set<ONDEXConcept> conceptview) {
		Map<String, Integer> map = new Hashtable<String, Integer>();
		String acc = null;
        int i = 0;
        int j = 0;
        int k = 0;
        for (ONDEXConcept c : conceptview) {
            i++;
	        Set<ConceptAccession> ci = c.getConceptAccessions();
	        for (ConceptAccession a : ci) {
	            j++;
	            CV conceptCV = a.getElementOf();
	            if (conceptCV.getId() == "UNIPROTKB" || conceptCV.getId() == "uniprot_id") {
	                k++;
                    acc = a.getAccession();
			        map.put(acc, c.getId());
	            }
	        }
		}
		return map;
	}
	private int getLastIdForConcept() {
        int stored = 0;
        Set<ONDEXConcept> view = graph.getConcepts();
        for (ONDEXConcept c : view) {
            Integer id = c.getId();
            if (id.intValue() > stored) {
                stored = id.intValue();
            }
        }
        return stored;
    }

    private Map<String, String> readAffyAnnotation(String filename) throws Exception {
    	Map<String, String> affyMapping = new Hashtable<String, String>();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line = null;
        int i = 0;
        while ((line = br.readLine()) != null) {
            String[] tokens = line.trim().split("\t");
            String[] affy_ids = tokens[1].split(":");
			for (String affy_id:affy_ids) {
				//this is doing lots of overwriting
                //use uniprot_id as key?
                String ups = null;
                try {
                    ups = affyMapping.get(affy_id);
                }
                catch (NullPointerException npe) {
                    affyMapping.put(affy_id, tokens[0]);
                }
                if (ups != null) {
                    ups = ups+":"+tokens[0];
                    affyMapping.put(affy_id, ups);
                }
                else {
                    affyMapping.put(affy_id, tokens[0]);
                }
			}
        }
        br.close();
        System.out.println(i+"******"+affyMapping.size()+"******");
		return affyMapping;
    }
    
	private Map<String, Float> readExpressionData(String filename, String tissue) throws Exception {
		Map<String, Float> expressionData = new Hashtable<String, Float>();
		int[] tissue_index = getColumnNumbers(filename, tissue);
		BufferedReader br = new BufferedReader(new FileReader(filename));
        String line = br.readLine(); //skip first line of the file (header)
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
			//Prostate expression vals in cols 93 & 94
			float total = 0.0f;
			for (int index:tissue_index) {
				Float value = Float.valueOf(tokens[index]);
				total = total + value.floatValue();
            }
			float expressionValue = total/(float)tissue_index.length;
			expressionData.put(tokens[0], new Float(expressionValue));
        }
		return expressionData;
    }

	private int[] getColumnNumbers(String f, String t) throws Exception {
		//assume 2 columns, because for default dataset there's always 2 values.
		int[] numbers = new int[2];
		BufferedReader br = new BufferedReader(new FileReader(f));
		String header = br.readLine();
		br.close();
		String[] tokens = header.split("\t");
		int i = 0;
		int j = 0;
		for (String token:tokens) {
			if (token.equalsIgnoreCase(t)) {
				numbers[j] = i;
				j++;
			}
			i++;
		}
		return numbers;
	}


    private String[] getExpressionInfo(String uniprot) {
        String affyIDs = affyMapping.get(uniprot);
		if (affyIDs != null) {
			return affyIDs.split(":");
		}
		else {
			return null;
		}
    }
}
		/*
		 * TISSUE TYPES:
		 * 1,2     ColorectalAdenocarcinoma
		 * 3,4     WHOLEBLOOD
		 * 5,6     BM-CD33+Myeloid
		 * 7,8     PB-CD14+Monocytes
		 * 9,10    PB-BDCA4+Dentritic_Cells
		 * 11,12   PB-CD56+NKCells
		 * 13,14   PB-CD4+Tcells
		 * 15,16   PB-CD8+Tcells
		 * 17,18   PB-CD19+Bcells
		 * 19,20   BM-CD105+Endothelial
		 * 21,22   BM-CD34+
		 * 23,24   leukemialymphoblastic(molt4)
		 * 25,26   721_B_lymphoblasts
		 * 27,28   lymphomaburkittsRaji
		 * 29,30   leukemiapromyelocytic(hl60)
		 * 31,32   lymphomaburkittsDaudi
		 * 33,34   leukemiachronicmyelogenous(k562)
		 * 35,36   thymus
		 * 37,38   Tonsil
		 * 39,40   lymphnode
		 * 41,42   fetalliver
		 * 43,44   BM-CD71+EarlyErythroid
		 * 45,46   bonemarrow
		 * 47,48   TemporalLobe
		 * 49,50   globuspallidus
		 * 51,52   CerebellumPeduncles
		 * 53,54   cerebellum
		 * 55,56   caudatenucleus
		 * 57,58   WholeBrain
		 * 59,60   ParietalLobe
		 * 61,62   MedullaOblongata
		 * 63,64   Amygdala
		 * 65,66   PrefrontalCortex
		 * 67,68   OccipitalLobe
		 * 69,70   Hypothalamus
		 * 71,72   Thalamus
		 * 73,74   subthalamicnucleus
		 * 75,76   CingulateCortex
		 * 77,78   Pons
		 * 79,80   spinalcord
		 * 81,82   fetalbrain
		 * 83,84   adrenalgland
		 * 85,86   Lung
		 * 87,88   Heart
		 * 89,90   Liver
		 * 91,92   kidney
		 * 93,94   Prostate
		 * 95,96   Uterus
		 * 97,98   Thyroid
		 * 99,100  fetalThyroid
		 * 101,102 fetallung
		 * 103,104 PLACENTA
		 * 105,106 CardiacMyocytes
		 * 107,108 SmoothMuscle
		 * 109,110 bronchialepithelialcells
		 * 111,112 ADIPOCYTE
		 * 113,114 Pancreas
		 * 115,116 PancreaticIslets
		 * 117,118 testis
		 * 119,120 TestisLeydigCell
		 * 121,122 TestisGermCell
		 * 123,124 TestisInterstitial
		 * 125,126 TestisSeminiferousTubule
		 * 127,128 salivarygland
		 * 129,130 trachea
		 * 131,132 AdrenalCortex
		 * 133,134 Ovary
		 * 135,136 Appendix
		 * 137,138 skin
		 * 139,140 ciliaryganglion
		 * 141,142 TrigeminalGanglion
		 * 143,144 atrioventricularnode
		 * 145,146 DRG
		 * 147,148 SuperiorCervicalGanglion
		 * 149,150 SkeletalMuscle
		 * 151,152 UterusCorpus
		 * 153,154 TONGUE
		 * 155,156 OlfactoryBulb
		 * 157,158 Pituitary
		 *
		 */
