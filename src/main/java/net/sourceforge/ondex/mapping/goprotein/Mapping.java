package net.sourceforge.ondex.mapping.goprotein;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.ArrayList;

import net.sourceforge.ondex.core.AttributeName;
import net.sourceforge.ondex.core.CV;
import net.sourceforge.ondex.core.ConceptAccession;
import net.sourceforge.ondex.core.ConceptClass;
import net.sourceforge.ondex.core.EvidenceType;
import net.sourceforge.ondex.core.GDS;
import net.sourceforge.ondex.core.ONDEXConcept;
import net.sourceforge.ondex.core.ONDEXGraph;
import net.sourceforge.ondex.core.ONDEXRelation;
import net.sourceforge.ondex.core.RelationType;
import net.sourceforge.ondex.mapping.ONDEXMapping;

import net.sourceforge.ondex.args.ArgumentDefinition;
import net.sourceforge.ondex.args.BooleanArgumentDefinition;
import net.sourceforge.ondex.args.StringArgumentDefinition;
import net.sourceforge.ondex.args.StringMappingPairArgumentDefinition;
//import net.sourceforge.ondex.mapping.MappingArguments;



public class Mapping extends ONDEXMapping {
    /**
     * Specifies neccessary arguments for this mapping.
     * 
     * @return ArgumentDefinition<?>[]
     */
    public ArgumentDefinition<?>[] getArgumentDefinitions() {
        return new ArgumentDefinition<?>[0];
    }
    public String getName() {
        return new String("GOA -> UniProt Mapping method");
    }
    public String getId() {
        return "goprotein";
    }
    public String getVersion() {
        return new String("30.07.2010");
    }
    public boolean requiresIndexedGraph() {
        return false;
    }
    public String[] requiresValidators() {
        return new String[0];
    }
    public void start() throws Exception{
        ArrayList goaConcepts = new ArrayList();
        ArrayList uniprotConcepts = new ArrayList();
        Set<ONDEXConcept> itConcepts = graph.getConcepts();
        for (ONDEXConcept concept: itConcepts) {
            CV conceptCV = concept.getElementOf();
            ConceptClass conceptc = concept.getOfType();
            String classId = conceptc.getId();
            if (classId == "Protein") {
            //add to kegg list
                if (conceptCV.getId() == "GOAEBI") {
                    goaConcepts.add(concept);
                }
            //add to uniprot list
                else if (conceptCV.getId().contains("UNIPROTKB")) {
                    uniprotConcepts.add(concept);
                }
            }
        }
        Object[] goaC = goaConcepts.toArray();
        Object[] uniprotC = uniprotConcepts.toArray();
        for (Object c: goaConcepts) {
			ONDEXConcept gc = (ONDEXConcept) c;
            Set<ConceptAccession> accessions = gc.getConceptAccessions();
            for (ConceptAccession accession: accessions) {
                String acc = accession.getAccession().trim();
                CV accCV = accession.getElementOf();
                if (accCV.getId() == "UNIPROTKB") {
                    //need to split this up if it is too long
                    if (acc.length() > 6) {
                        int i = 0;
                        int j = 6;
                        while (j <= acc.length()) {
                            String subacc = acc.substring(i, j);
                            System.out.println(subacc);
                            matchAccessions(uniprotC, subacc, gc);
                            i += 6;
                            j += 6;
                        }
                    }
                    else {
                        matchAccessions(uniprotC, acc, gc);
                    }
                }
            }
        }
    }
    public void matchAccessions(Object[] uniprotArray, String accession, ONDEXConcept toConcept) throws Exception {
        RelationType equRelation = requireRelationType("equ");
        EvidenceType accET = requireEvidenceType("ACC");
		for (Object uc: uniprotArray) {
			ONDEXConcept upc = (ONDEXConcept) uc;
            Set<ConceptAccession> up_accessions = upc.getConceptAccessions();
            for (ConceptAccession up_accession : up_accessions) {
                String u_a = up_accession.getAccession().trim();
                CV uaccCV = up_accession.getElementOf();
                if (uaccCV.getId() == "UNIPROTKB") {
                    if (accession.matches(u_a)) {
                        ONDEXRelation r1 = graph.getFactory().createRelation(upc, toConcept, equRelation, accET);
                    }
                }
            }
        }
    }
}
