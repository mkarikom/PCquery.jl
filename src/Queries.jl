# strings for sample queries

function path_pc(pathVals::String)
    """
    PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
    PREFIX dcterms: <http://purl.org/dc/terms/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>

    SELECT DISTINCT ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller ?transcriptionController ?cxRefL ?cxRefR ?upRefL ?upRefR ?upDnaRefL ?upDnaRefR
    # (GROUP_CONCAT(?leref; SEPARATOR = ",") AS ?lerefs) (GROUP_CONCAT(?reref; SEPARATOR = ",") AS ?rerefs)

    WHERE {
        OPTIONAL { ?left a bp:Complex ;
                         bp:xref [ a bp:UnificationXref ;
                                   bp:db "reactome"^^xsd:string ;
                                   bp:id ?cxRefL] }

        OPTIONAL { ?right a bp:Complex ;
                         bp:xref [ a bp:UnificationXref ;
                                   bp:db "reactome"^^xsd:string ;
                                   bp:id ?cxRefR] }

        OPTIONAL { ?left a bp:Protein ;
                         bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
                                                                               bp:db ?dbname ;
                                                                               bp:id ?upRefL ]
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Protein ;
                          bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
                                                                                bp:db ?dbname ;
                                                                                bp:id ?upRefR ]
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?left a bp:Protein ;
                         bp:entityReference [ a bp:ProteinReference ;
                                              bp:xref [ bp:db ?dbname ;
                                                        bp:id ?upRefL ] ] .
                    FILTER ( NOT EXISTS {?left bp:memberPhysicalEntity ?lmem } )
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Protein ;
                          bp:entityReference [ a bp:ProteinReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:id ?upRefR ] ] .
                    FILTER ( NOT EXISTS {?right bp:memberPhysicalEntity ?rmem } )
                    FILTER regex(?dbname, "^.*uniprot") }

        # single genes
        OPTIONAL { ?left a bp:Dna ;
                          bp:entityReference [ a bp:DnaReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:id ?upDnaRefL] ]
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Dna ;
                          bp:entityReference [ a bp:DnaReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:id ?upDnaRefR] ]
                    FILTER regex(?dbname, "^.*uniprot") }

        # gene family
        OPTIONAL { ?left a bp:Dna ;
                          bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:RelationshipXref ;
                                                                                bp:db ?dbname ;
                                                                                bp:id ?upDnaRefL]
                          FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Dna ;
                          bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:RelationshipXref ;
                                                                                bp:db ?dbname ;
                                                                                bp:id ?upDnaRefR]
                          FILTER regex(?dbname, "^.*uniprot") }
        {
            SELECT ?pw ?pathname ?pc ?ctype
                    ?left ?right ?ltype ?rtype ?llocref ?rlocref
                    ?cltype ?controller ?transcriptionController
            WHERE {
                {# the actual values of left and right (we will not recognize "left side DNA" in downstream calculations)
                    ?pw bp:name ?pathname ;
                        bp:pathwayComponent* ?pc .
                     # VALUES ?pw {<http://identifiers.org/reactome/R-HSA-2122948>
                     #             <http://identifiers.org/reactome/R-HSA-2122947>}
                    VALUES ?pw $pathVals

                    ?pc a ?ctype .
                    ?pc a bp:Interaction .
                    ?left ^bp:left ?pc ;
                           a ?ltype ;
                           bp:cellularLocation [ bp:term ?lloc ;
                                                 bp:xref ?llocref ] ;
                           bp:displayName ?lname .
                    ?right ^bp:right ?pc ;
                           a ?rtype ;
                           bp:cellularLocation [ bp:term ?rloc ;
                                                 bp:xref ?rlocref ] ;
                           bp:displayName ?rname .
                    OPTIONAL { ?pc ^bp:controlled [ bp:controlType ?cltype ;
                                                    bp:controller ?controller ]}
                    VALUES ?rtype {bp:Protein bp:Complex}
                    VALUES ?ltype {bp:Protein bp:Complex} }
                UNION
                {# we count reactions where left2 is dna, left2 is an activator complex, and right is the trans init complex, then report a new reaction where left is the trans init complex and right is the protein encoded by the gene
                    ?pw bp:name ?pathname ;
                        bp:pathwayComponent* ?pc .
                     # VALUES ?pw {<http://identifiers.org/reactome/R-HSA-2122948>
                     #             <http://identifiers.org/reactome/R-HSA-2122947>}

                    VALUES ?pw $pathVals

                    ?pc a ?ctype .
                    ?pc a bp:BiochemicalReaction ; # identify the controlling complex for a transcriptional event
                        bp:left ?right ;
                        ^bp:controlled [ a bp:Control ;
                                         bp:controlType ?cltype ;
                                         bp:controller ?left ;
                                         bp:xref [ a bp:UnificationXref ;
                                                   bp:db ?dbname ;
                                                   bp:id ?transcriptionController ] ] .
                    ?left a ?ltype ;
                           bp:cellularLocation [ bp:term ?lloc ;
                                                 bp:xref ?llocref ] .

                    ?right a bp:Dna ;
                           bp:cellularLocation [ bp:term ?rloc ;
                                                 bp:xref ?rlocref ] .

                    VALUES ?rtype {bp:Protein}
                    VALUES ?ltype {bp:Protein bp:Complex} }
                    VALUES ?ctype {bp:BiochemicalReaction} } }
    }

    GROUP BY ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller ?transcriptionController ?cxRefL ?cxRefR ?upRefL ?upRefR ?upDnaRefL ?upDnaRefR
    """
end

function path_pc_recurse(pathVals::String)
    """
    PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
    PREFIX dcterms: <http://purl.org/dc/terms/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>

    SELECT DISTINCT ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller ?transcriptionController ?cxRefL ?cxRefR ?upRefL ?upRefR ?upDnaRefL ?upDnaRefR
    (GROUP_CONCAT(?cxPrimaryProtL; SEPARATOR = ",") AS ?cxPrimaryProtLcat)
    (GROUP_CONCAT(?cxPrimaryProtR; SEPARATOR = ",") AS ?cxPrimaryProtRcat)
    (GROUP_CONCAT(?cxPrimaryProtLocL; SEPARATOR = ",") AS ?cxPrimaryProtLocLcat)
    (GROUP_CONCAT(?cxPrimaryProtLocR; SEPARATOR = ",") AS ?cxPrimaryProtLocRcat)
    (GROUP_CONCAT(?cxPrimaryCxL; SEPARATOR = ",") AS ?cxPrimaryCxLcat)
    (GROUP_CONCAT(?cxPrimaryCxR; SEPARATOR = ",") AS ?cxPrimaryCxRcat)
    (GROUP_CONCAT(?cxPrimaryCxLocL; SEPARATOR = ",") AS ?cxPrimaryCxLocLcat)
    (GROUP_CONCAT(?cxPrimaryCxLocR; SEPARATOR = ",") AS ?cxPrimaryCxLocRcat)

    WHERE {
    		## complexes: primary protein components
    		OPTIONAL { ?left a bp:Complex ;
    										 bp:xref [ a bp:UnificationXref ;
    															 bp:db "reactome"^^xsd:string ;
    															 bp:id ?cxRefL] ;
    															 bp:component [ a bp:Protein ;
    																							bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
    																																																		bp:db ?dbname ;
    																																																		bp:id ?cxPrimaryProtL] ;
    																							bp:cellularLocation [ bp:xref ?cxPrimaryProtLocL ] ] .
    							FILTER( regex(?dbname, "^.*uniprot") && BOUND(?cxRefL) && BOUND(?cxPrimaryProtL) && BOUND(?cxPrimaryProtLocL) )}

    		OPTIONAL { ?right a bp:Complex ;
    										 bp:xref [ a bp:UnificationXref ;
    															 bp:db "reactome"^^xsd:string ;
    															 bp:id ?cxRefR] ;
    															 bp:component [ a bp:Protein ;
    																							bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
    																																																		bp:db ?dbname ;
    																																																		bp:id ?cxPrimaryProtR] ;
    																							bp:cellularLocation [ bp:xref ?cxPrimaryProtLocR ] ] .
    							FILTER( regex(?dbname, "^.*uniprot") && BOUND(?cxRefR) && BOUND(?cxPrimaryProtR) && BOUND(?cxPrimaryProtLocR) )  }

    		 OPTIONAL { ?left a bp:Complex ;
    										 bp:xref [ a bp:UnificationXref ;
    															 bp:db "reactome"^^xsd:string ;
    															 bp:id ?cxRefL] ;
    															 bp:component [ a bp:Protein ;
    																							bp:entityReference/bp:xref [ a bp:UnificationXref ;
    																																					 bp:db ?dbname ;
    																																					 bp:id ?cxPrimaryProtL] ;
    																							bp:cellularLocation [ bp:xref ?cxPrimaryProtLocL ] ] .
    							 FILTER( regex(?dbname, "^.*uniprot") && BOUND(?cxRefL) && BOUND(?cxPrimaryProtL) && BOUND(?cxPrimaryProtLocL) )  }

    		 OPTIONAL { ?right a bp:Complex ;
    										 bp:xref [ a bp:UnificationXref ;
    															 bp:db "reactome"^^xsd:string ;
    															 bp:id ?cxRefR] ;
    															 bp:component [ a bp:Protein ;
    																							bp:entityReference/bp:xref [ a bp:UnificationXref ;
    																																					 bp:db ?dbname ;
    																																					 bp:id ?cxPrimaryProtR] ;
    																						  bp:cellularLocation [ bp:xref ?cxPrimaryProtLocR ] ] .
    							 FILTER( regex(?dbname, "^.*uniprot") && BOUND(?cxRefR) && BOUND(?cxPrimaryProtR) && BOUND(?cxPrimaryProtLocR) )  }

    		 # complexes: primary complex components
    		 OPTIONAL { ?left a bp:Complex ;
     										 bp:xref [ a bp:UnificationXref ;
     															 bp:db "reactome"^^xsd:string ;
     															 bp:id ?cxRefL] ;
     															 bp:component [ a bp:Complex ;
     																							bp:xref [ a bp:UnificationXref ;
    																							bp:db "reactome"^^xsd:string ;
    																							bp:id ?cxPrimaryCxL] ;
     																							bp:cellularLocation [ bp:xref ?cxPrimaryCxLocL ] ] .
     							FILTER( regex(?dbname, "^.*uniprot") && BOUND(?cxRefL) && BOUND(?cxPrimaryCxL) && BOUND(?cxPrimaryCxLocL) )}

    		OPTIONAL { ?right a bp:Complex ;
    										 bp:xref [ a bp:UnificationXref ;
    															 bp:db "reactome"^^xsd:string ;
    															 bp:id ?cxRefR] ;
    															 bp:component [ a bp:Complex ;
    																							bp:xref [ a bp:UnificationXref ;
    																							bp:db "reactome"^^xsd:string ;
    																							bp:id ?cxPrimaryCxR] ;
    																							bp:cellularLocation [ bp:xref ?cxPrimaryCxLocR ] ] .
    							FILTER( regex(?dbname, "^.*uniprot") && BOUND(?cxRefR) && BOUND(?cxPrimaryCxR) && BOUND(?cxPrimaryCxLocR) )}

    		# proteins
    		OPTIONAL { ?left a bp:Protein ;
    										 bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
    																																					 bp:db ?dbname ;
    																																					 bp:id ?upRefL ]
    								FILTER regex(?dbname, "^.*uniprot") }

    		OPTIONAL { ?right a bp:Protein ;
    											bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
    																																						bp:db ?dbname ;
    																																						bp:id ?upRefR ]
    								FILTER regex(?dbname, "^.*uniprot") }

    		OPTIONAL { ?left a bp:Protein ;
    										 bp:entityReference [ a bp:ProteinReference ;
    																					bp:xref [ bp:db ?dbname ;
    																										bp:id ?upRefL ] ] .
    								FILTER ( NOT EXISTS {?left bp:memberPhysicalEntity ?lmem } )
    								FILTER regex(?dbname, "^.*uniprot") }

    		OPTIONAL { ?right a bp:Protein ;
    											bp:entityReference [ a bp:ProteinReference ;
    																					 bp:xref [ bp:db ?dbname ;
    																										 bp:id ?upRefR ] ] .
    								FILTER ( NOT EXISTS {?right bp:memberPhysicalEntity ?rmem } )
    								FILTER regex(?dbname, "^.*uniprot") }

    		# dna: single genes
    		OPTIONAL { ?left a bp:Dna ;
    											bp:entityReference [ a bp:DnaReference ;
    																					 bp:xref [ bp:db ?dbname ;
    																										 bp:id ?upDnaRefL] ]
    								FILTER regex(?dbname, "^.*uniprot") }

    		OPTIONAL { ?right a bp:Dna ;
    											bp:entityReference [ a bp:DnaReference ;
    																					 bp:xref [ bp:db ?dbname ;
    																										 bp:id ?upDnaRefR] ]
    								FILTER regex(?dbname, "^.*uniprot") }

    		# dna: gene family
    		OPTIONAL { ?left a bp:Dna ;
    											bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:RelationshipXref ;
    																																						bp:db ?dbname ;
    																																						bp:id ?upDnaRefL]
    											FILTER regex(?dbname, "^.*uniprot") }

    		OPTIONAL { ?right a bp:Dna ;
    											bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:RelationshipXref ;
    																																						bp:db ?dbname ;
    																																						bp:id ?upDnaRefR]
    											FILTER regex(?dbname, "^.*uniprot") }
    		{
    				SELECT ?pw ?pathname ?pc ?ctype
    								?left ?right ?ltype ?rtype ?llocref ?rlocref
    								?cltype ?controller ?transcriptionController
    				WHERE {
    						{# the actual values of left and right (we will not recognize "left side DNA" in downstream calculations)
    								?pw bp:name ?pathname ;
    										bp:pathwayComponent* ?pc .
    									VALUES ?pw {<http://identifiers.org/reactome/R-HSA-2122948>
    															<http://identifiers.org/reactome/R-HSA-2122947>}
    #                    VALUES ?pw $pathVals

    								?pc a ?ctype .
    								?pc a bp:Interaction .
    								?left ^bp:left ?pc ;
    											 a ?ltype ;
    											 bp:cellularLocation [ bp:term ?lloc ;
    																						 bp:xref ?llocref ] ;
    											 bp:displayName ?lname .
    								?right ^bp:right ?pc ;
    											 a ?rtype ;
    											 bp:cellularLocation [ bp:term ?rloc ;
    																						 bp:xref ?rlocref ] ;
    											 bp:displayName ?rname .
    								OPTIONAL { ?pc ^bp:controlled [ bp:controlType ?cltype ;
    																								bp:controller ?controller ]}
    								VALUES ?rtype {bp:Protein bp:Complex}
    								VALUES ?ltype {bp:Protein bp:Complex}
    #                    VALUES ?rtype {bp:Complex}
    #                    VALUES ?ltype {bp:Complex}
    				}
    						UNION
    						{# count reactions where left2 is dna, left2 is an activator complex, and right is the trans init complex, then report a new reaction where left is the trans init complex and right is the protein encoded by the gene
    								?pw bp:name ?pathname ;
    										bp:pathwayComponent* ?pc .
    									VALUES ?pw {<http://identifiers.org/reactome/R-HSA-2122948>
    															<http://identifiers.org/reactome/R-HSA-2122947>}

    #                    VALUES ?pw $pathVals

    								?pc a ?ctype .
    								?pc a bp:BiochemicalReaction ; # identify the controlling complex for a transcriptional event
    										bp:left ?right ;
    										^bp:controlled [ a bp:Control ;
    																		 bp:controlType ?cltype ;
    																		 bp:controller ?left ;
    																		 bp:xref [ a bp:UnificationXref ;
    																							 bp:db ?dbname ;
    																							 bp:id ?transcriptionController ] ] .
    								?left a ?ltype ;
    											 bp:cellularLocation [ bp:term ?lloc ;
    																						 bp:xref ?llocref ] .

    								?right a bp:Dna ;
    											 bp:cellularLocation [ bp:term ?rloc ;
    																						 bp:xref ?rlocref ] .

    								VALUES ?rtype {bp:Protein}
    								VALUES ?ltype {bp:Protein bp:Complex} }
    								VALUES ?ctype {bp:BiochemicalReaction} } }

    #    FILTER ( !BOUND( ?cxPrimaryProtR ) )
    }

    GROUP BY ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller ?transcriptionController ?cxRefL ?cxRefR ?upRefL ?upRefR ?upDnaRefL ?upDnaRefR
    """
end

function np_annotations(goFilter::String,
                        upFilter::String)
    """
    PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
    PREFIX dcterms: <http://purl.org/dc/terms/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
    PREFIX up: <http://purl.uniprot.org/uniprot/>
    PREFIX np: <http://nextprot.org/rdf#>
    PREFIX cv: <http://nextprot.org/rdf/terminology/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX unipage: <http://www.uniprot.org/uniprot/>

    SELECT DISTINCT *
    WHERE {
            ?entry a np:Entry ;
                   np:swissprotPage ?uniprot ;
                   np:isoform ?isoform .
            ?isoform np:goMolecularFunction/np:term/np:childOf ?goterm ;
                     np:goMolecularFunction/np:term ?term
            VALUES ?goterm $goFilter

            VALUES ?uniprot $upFilter
    } LIMIT 100
    """
end
