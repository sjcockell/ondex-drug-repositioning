package net.sourceforge.ondex.parser.ethera_data;

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
import net.sourceforge.ondex.core.ConceptAccession;
import net.sourceforge.ondex.core.ConceptName;

import net.sourceforge.ondex.exception.type.MetaDataMissingException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Writer;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Hashtable;
import java.util.Map;
import java.util.ArrayList;
import java.util.Iterator;

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
	private CV uniprotCV;
    private CV hprdCV;
    private CV hugoCV;
    private CV affyCV;
    
    private CV newCV;
    private ConceptClass newCC;
    private EvidenceType newET;
    private RelationType newRelation;
    
    private ONDEXConcept MatchedConcept1 = null;
    private ONDEXConcept MatchedConcept2 = null;
    private ArrayList matchedConcepts = new ArrayList();
   
    public static final String INPUT_FILE_DESC = "File to import";
    private ArgumentDefinition<?>[] argdefs = new ArgumentDefinition<?>[]{
        new FileArgumentDefinition(FileArgumentDefinition.INPUT_FILE,
            INPUT_FILE_DESC, true, true, false, false)};

    public boolean readsDirectory() {
        return false;
    }

    public boolean readsFile() {
        return true;
    }

    public ArgumentDefinition<?>[] getArgumentDefinitions() {
        /*
        new ArgumentDefinition<?>[] {
        new StringArgumentDefinition(params)
        }
         */
        return argdefs;
    }

    public String getName() {
        return "ETherapeutics Data Parser";
    }

    public String getId() {
        return "ethera_data";
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
        ArrayList proteins = getProteinsFromNetwork();
        //get file name
        String filename = (String)pa.getUniqueValue(FileArgumentDefinition.INPUT_FILE);
        BufferedReader br = new BufferedReader(new FileReader(filename));
        //process the file, line at a time
        //output files
        Writer output1 = new BufferedWriter(new FileWriter("unmatched.txt"));
        //Writer output2 = new BufferedWriter(new FileWriter(outfile2));
        String line = null;
        int i = 0;
        while ((line = br.readLine()) != null) {
            ONDEXRelation r1 = null;
            String[] columns = line.split(",");
            //Setup concept one (<from> protein)
            /*    
            METHODS FOR ANNOTATING CONCEPTS
            createConceptAccession(String accession, CV, boolean ambiguous)
            createConceptName(String name, boolean preferred
            createConceptGDS(AttributeName attrname, Object value, boolean doIndex)

             */
            String cID1 = columns[0].trim();
            String cID2 = columns[1].trim();
            MatchedConcept1 = null;
            MatchedConcept2 = null;
            for (int index = 0; index < proteins.size(); index += 1) {
                ONDEXConcept protein = (ONDEXConcept) proteins.get(index);
                //this matches for HPRD accessions in the file
                boolean match1 = checkEquivalence(cID1, protein, i, 1);
                if (match1 == true) {
                    i++;
                }
                boolean match2 = checkEquivalence(cID2, protein, i, 2);
                if (match2 == true) {
                    i++;
                }
            }
            if (MatchedConcept1 != null && MatchedConcept2 != null) {
                r1 = graph.getFactory().createRelation(MatchedConcept1, MatchedConcept2, proteinProteinRelation, newET);
            }
            if (MatchedConcept1 == null) {
                output1.write(cID1+"\n");
            }
            if (MatchedConcept2 == null) {
                output1.write(cID2+"\n");
            }
        }
        System.out.println(i+" matches made");
    }

    private ArrayList getProteinsFromNetwork() {
        ArrayList uniprotConcepts = new ArrayList();
        Iterator<ONDEXConcept> itConcepts = graph.getConcepts().iterator();
        while (itConcepts.hasNext()) {
            ONDEXConcept concept = itConcepts.next();
            CV conceptCV = concept.getElementOf();
            ConceptClass conceptc = concept.getOfType();
            String classId = conceptc.getId();
            if (classId == "Protein") {
                if (conceptCV.getId().contains("UNIPROTKB")) {
                    uniprotConcepts.add(concept);
                }
            }
        }
        return uniprotConcepts;
    }

    private boolean checkEquivalence(String id, ONDEXConcept pro, int count, int which) {
        ONDEXConcept c1 = null;
        ONDEXRelation r1 = null;
        boolean match = false;
        if (id.length() > 0) {
        if (id.startsWith("HPRD")) {
            //simple accession mapping
            String accession1 = id.substring(5);
            ConceptAccession acc1 = pro.getConceptAccession(accession1, hprdCV);
            if (acc1 != null) {
                match = true;
                c1 = checkExists(acc1.getAccession(), 1);
                //make new concept, link to the protein above
                if (c1 == null) {
                    c1 = graph.getFactory().createConcept(++count+"", newCV, newCC, newET);
                    c1.createConceptName(id, false);
                    c1.createConceptAccession(id, newCV, false);
                    matchedConcepts.add(c1);
                }
                r1 = graph.getFactory().createRelation(c1, pro, newRelation, newET);
            }
        }
        else {
            //String matching - CASE INSENSITIVE??
            ConceptName name1 = pro.getConceptName(id);
            if (name1 != null) {
                c1 = checkExists(name1.getName(), 2);
                match = true;
                if (c1 == null) {
                c1 = graph.getFactory().createConcept(++count+"", newCV, newCC, newET);
                c1.createConceptName(id, false);
                matchedConcepts.add(c1);
                }
                r1 = graph.getFactory().createRelation(c1, pro, newRelation, newET);
            }

        }
        }
        if (which == 1) {
            MatchedConcept1 = c1;
        }
        else if (which == 2) {
            MatchedConcept2 = c1;
        }
        return match;
    }

    private ONDEXConcept checkExists(String id, int type) {
        ONDEXConcept toReturn = null;
        boolean matches = false;
        if (type == 1) {
            for (int index = 0; index < matchedConcepts.size(); index += 1) {
                ONDEXConcept c = (ONDEXConcept) matchedConcepts.get(index);
                ConceptAccession acc = c.getConceptAccession(id, newCV);
                if (acc != null) {
                    matches = true;
                    toReturn = c;
                }
            }
        }
        else if (type == 2) {
            for (int index = 0; index < matchedConcepts.size(); index += 1) {
                ONDEXConcept c = (ONDEXConcept) matchedConcepts.get(index);
                ConceptName acc = c.getConceptName(id);
                if (acc != null) {
                    matches = true;
                    toReturn = c;
                }
            }
        }
        return toReturn;
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

        newCV = requireCV("eTherapeutics");
        newCC = requireConceptClass("eTherapeutics_Protein");
        newET = requireEvidenceType("ETANN");
        newRelation = requireRelationType("equ");
    }

}


