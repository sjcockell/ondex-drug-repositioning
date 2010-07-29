package net.sourceforge.ondex.mapping.drugmap;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.ArrayList;

import net.sourceforge.ondex.mapping.ONDEXMapping;

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
        return new String("DrugBank Mapping method");
    }
    public String getId() {
        return "kegguni";
    }
    public String getVersion() {
        return new String("02.12.2009");
    }
    public boolean requiresIndexedGraph() {
        return false;
    }
    public String[] requiresValidators() {
        return new String[0];
    }
    public void start() throws Exception{
        ArrayList drugConcepts = new ArrayList();
        ArrayList uncConcepts = new ArrayList();
        Set<ONDEXConcept> itConcepts = graph.getConcepts();
        for (ONDEXConcept concept: itConcepts) {
            CV conceptCV = concept.getElementOf();
            ConceptClass conceptc = concept.getOfType();
            String classId = conceptc.getId();
            if (classId == "Comp") {
            //add to kegg list
                    drugConcepts.add(concept);
            }
            //add to uniprot list
            else if (classId.equals("Compound")) {
                    uncConcepts.add(concept);
            }
        }
        Object[] drugC = drugConcepts.toArray();
        Object[] uncC = uncConcepts.toArray();
        for (Object c: uncC) {
            ONDEXConcept UncConcept = (ONDEXConcept) c;
            Set<ConceptAccession> accessions = UncConcept.getConceptAccessions();
            for (ConceptAccession accession : accessions) {
                String acc = accession.getAccession().trim();
                CV accCV = accession.getElementOf();
                if (accCV.getId().equals("drugbank_id")) {
                    //need to split this up if it is too long
                    System.out.println(acc);
                    matchAccessions(drugC, acc, UncConcept);
                }
            }
        }
    }
    public void matchAccessions(Object[] drugArray, String accession, ONDEXConcept toConcept) throws Exception {
        RelationType equRelation = requireRelationType("temp");
        EvidenceType accET = requireEvidenceType("ACC");
        for (Object dc: drugArray) {
            ONDEXConcept DrugConcept = (ONDEXConcept) dc;
            Set<ConceptAccession> dr_accessions = DrugConcept.getConceptAccessions();
            for (ConceptAccession dr_accession: dr_accessions) {
                String d_a = dr_accession.getAccession().trim();
                CV daccCV = dr_accession.getElementOf();
                if (daccCV.getId() == "DRUGBANK") {
                    if (accession.matches(d_a)) {
                        ONDEXRelation r1 = graph.getFactory().createRelation(DrugConcept, toConcept, equRelation, accET);
                    }
                }
            }
        }
    }
}
