package net.sourceforge.ondex.parser.etherapeutics;

import net.sourceforge.ondex.parser.ONDEXParser;
import net.sourceforge.ondex.args.ArgumentDefinition;
import net.sourceforge.ondex.args.FileArgumentDefinition;
import net.sourceforge.ondex.core.ONDEXConcept;
import net.sourceforge.ondex.core.ONDEXRelation;
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

public class Parser extends ONDEXParser {
    /* 
    Author - sjcockell
    Date - 20090630
    Version - 0.0.1
     */

    private ConceptClass proteinCC;
    private RelationType proteinProteinRelation;
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
	private Map<String, Float> expressionData = new Hashtable<String, Float>();
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
            new FileArgumentDefinition(FileArgumentDefinition.INPUT_DIR, "Directory with input files", true, true, true, false)
        };
    }

    public String getName() {
        return "ETherapeutics";
    }

    public String getId() {
        return "etherapeutics";
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
        //get directory name
        String dirname = (String)pa.getUniqueValue(FileArgumentDefinition.INPUT_DIR);
        processProteins(dirname + "/proteins");
        readAffyAnnotation(dirname + "/uniprot_affy");
		readExpressionData(dirname + "/expression");
        //read primary file
        String filename = dirname + "/interactions";
        BufferedReader br = new BufferedReader(new FileReader(filename));
        //process the file, line at a time
        String line = null;
        int i = 0;
        while ((line = br.readLine()) != null) {
            ONDEXConcept c1 = null;
            ONDEXConcept c2 = null;
            ONDEXRelation r1 = null;
            String[] columns = line.split("\t");
            //Setup concept one (<from> protein)
            /*    
            METHODS FOR ANNOTATING CONCEPTS
            createConceptAccession(String accession, CV, boolean ambiguous)
            createConceptName(String name, boolean preferred
            createConceptGDS(AttributeName attrname, Object value, boolean doIndex)

             */
            String cID1 = columns[0].trim();
            String cID2 = columns[1].trim();
            c1 = proteins.get(cID1);
            if (c1 == null) {
                c1 = graph.getFactory().createConcept(++i + "", genericCV, proteinCC, genericET);
                c1.createConceptName(cID1, false);
                c1.createConceptAccession(cID1, genericCV, false);
                String up = uniprotInfo.get(cID1);
                if (up != null) {
                    c1.createConceptAccession(up, uniprotCV, false);
                    String[] affyIDList = getExpressionInfo(up);
					if (affyIDList != null)  {
						float avExpression = calcAverageExpression(affyIDList);
						for (int trace = 0;trace < affyIDList.length;trace++) {
							c1.createConceptAccession(affyIDList[trace], affyCV, false);
						}
						c1.createGDS(expressionGDS, new Float(avExpression), false);
					}
                }
                String hprd = hprdInfo.get(cID1);
                if (hprd != null) {
                    c1.createConceptAccession(hprd, hprdCV, false);
                }
                String hugo = hugoInfo.get(cID1);
                if (hugo != null) {
                    c1.createConceptAccession(hugo, hugoCV, false);
                }
                String genesymbol = gsInfo.get(cID1);
                if (genesymbol != null) {
                    c1.createConceptName(genesymbol, false);
                }
                String genename = gnInfo.get(cID1);
                if (genename != null) {
                    c1.createConceptName(genename, true);
                }
                proteins.put(cID1, c1);
            }
            //Setup concept two (<to> protein)
            c2 = proteins.get(cID2);
            if (c2 == null) {
                c2 = graph.getFactory().createConcept(++i + "", genericCV, proteinCC, genericET);
                c2.createConceptAccession(cID2, genericCV, false);
                c2.createConceptName(cID2, true);
                String up = uniprotInfo.get(cID2);
                if (up != null) {
                    c2.createConceptAccession(up, uniprotCV, false);
                    String[] affyIDList = getExpressionInfo(up);
					if (affyIDList != null)  {
						float avExpression = calcAverageExpression(affyIDList);
						for (int trace = 0;trace < affyIDList.length;trace++) {
							c2.createConceptAccession(affyIDList[trace], affyCV, false);
						}
						c2.createGDS(expressionGDS, new Float(avExpression), false);
					}
                }
                String hprd = hprdInfo.get(cID2);
                if (hprd != null) {
                    c2.createConceptAccession(hprd, hprdCV, false);
                }
                String hugo = hugoInfo.get(cID2);
                if (hugo != null) {
                    c2.createConceptAccession(hugo, hugoCV, false);
                }
                String genesymbol = gsInfo.get(cID2);
                if (genesymbol != null) {
                    c2.createConceptName(genesymbol, false);
                }
                String genename = gnInfo.get(cID2);
                if (genename != null) {
                    c2.createConceptName(genename, true);
                }
                proteins.put(cID2, c2);
            }
            //Add relation, and annotate with G-Sesame Score
            r1 = graph.getFactory().createRelation(c2, c1, proteinProteinRelation, genericET);
            if (columns[2].trim().matches("None")) {
            } else {
                r1.createGDS(scoreGDS, new Float(columns[2].trim()), false);
            }
        }
    }

    private void fetchMetadata() throws MetaDataMissingException {
        proteinCC = requireConceptClass("Protein");
        proteinProteinRelation = requireRelationType("interacts_with");
        genericCV = requireCV("eTherapeutics");
        uniprotCV = requireCV("UNIPROTKB");
        hprdCV = requireCV("HPRD");
        hugoCV = requireCV("HGNC");
		affyCV = requireCV("AFFYMETRIX");
        genericET = requireEvidenceType("IMPD");
        scoreGDS = requireAttributeName("G-Sesame Score");
		expressionGDS = requireAttributeName("Average Expression");
    }

    private void processProteins(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        //process the file, line at a time
        String line = null;
        while ((line = br.readLine()) != null) {
            ProteinInfo pi = new ProteinInfo(line);
            if (pi.getUPID() != null) {
                uniprotInfo.put(pi.getETID(), pi.getUPID());
            }
            if (pi.getHPRDID() != null) {
                hprdInfo.put(pi.getETID(), pi.getHPRDID());
            }
            if (pi.getHUGOID() != null) {
                hugoInfo.put(pi.getETID(), pi.getHUGOID());
            }
            if (pi.getGeneSymbol() != null) {
                gsInfo.put(pi.getETID(), pi.getGeneSymbol());
            }
            if (pi.getGeneName() != null) {
                gnInfo.put(pi.getETID(), pi.getGeneName());
            }
        }
    }

    private void readAffyAnnotation(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line = null;
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            affyMapping.put(tokens[0], tokens[1]);
        }
    }
    
	private void readExpressionData(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line = br.readLine(); //skip first line of the file (header)
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
			//Prostate expression vals in cols 93 & 94
			Float value1 = Float.valueOf(tokens[93]);
			Float value2 = Float.valueOf(tokens[94]);
            float expressionValue = (value1.floatValue() + value2.floatValue())/2.0f;
			expressionData.put(tokens[0], new Float(expressionValue));
        }
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

	private float calcAverageExpression(String[] list) {
		Float total = new Float(0);
		for (int i = 0; i < list.length; i++) {
			Float exp = expressionData.get(list[0]);
			total = total + exp;
		}
		float av = total.floatValue() / (float)list.length;
		return av;
	}
}

class ProteinInfo {

    public String etID = null;
    public String upID = null;
    public String hprdID = null;
    public ArrayList goCCIDs = null;
    public ArrayList goMFIDs = null;
    public ArrayList goBPIDs = null;
    public String hugoID = null;
    public String geneSymbol = null;
    public String geneName = null;

    public ProteinInfo(String input) {
        String[] fields = new String[13];
        String[] info = input.split("\t");
        for (int i = 0; i < 13; i++) {
            try {
                if (info[i].equals("")) {
                    fields[i] = null;
                } else {
                    fields[i] = info[i];
                }
            } catch (ArrayIndexOutOfBoundsException e) {
                fields[i] = null;
            }
        }
        try {
            etID = fields[0].trim();
            upID = fields[10].trim();
            hprdID = fields[3].trim();
            hugoID = fields[4].trim();
            geneSymbol = fields[5].trim();
            geneName = fields[7].trim();
        } catch (NullPointerException ne) {
        }
    }

    public String getUPID() {
        return upID;
    }

    public String getETID() {
        return etID;
    }

    public String getHPRDID() {
        return hprdID;
    }

    public String getHUGOID() {
        return hugoID;
    }

    public String getGeneSymbol() {
        return geneSymbol;
    }

    public String getGeneName() {
        return geneName;
    }
}

