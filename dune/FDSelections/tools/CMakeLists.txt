simple_plugin(LongestRecoTrackSelector "tool"
    	  	        lardataobj_RecoBase
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas
 	  	           art_Utilities canvas

             )


simple_plugin(LongestRecoVertexTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas
             )

simple_plugin(HighestCTPMuonScoreRecoVertexTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas
    
                dune_TrackPID_tf
                dune_TrackPID_algorithms
                dune_TrackPID_products
             )

simple_plugin(HighestPandizzleScoreRecoVertexTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandizzleAlg    
                #dune_FDSelections_pandizzle
             )

simple_plugin(HighestPandrizzleScoreRecoVertexShowerSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandrizzleAlg    
                #dune_FDSelections_pandizzle
             )

simple_plugin(HighestPandrizzleScoreRecoVertexShowerSelectorCheatCut "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   FDSelection
                   PandrizzleAlg

                #dune_FDSelections_pandizzle
             )

simple_plugin(CheatRecoVertexTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandizzleAlg    
                #dune_FDSelections_pandizzle
             )

simple_plugin(CheatRecoTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandizzleAlg    
                #dune_FDSelections_pandizzle
             )

simple_plugin(CheatRecoPfoTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandizzleAlg    
                #dune_FDSelections_pandizzle
             )

simple_plugin(CheatRecoVertexShowerSelector "tool"
	  	        larcore_Geometry_Geometry_service
		          larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
 	  	       larsim_MCCheater_BackTrackerService_service
         	  larsim_MCCheater_ParticleInventoryService_service
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

            larsim_Utils
                   PandrizzleAlg    
             )

simple_plugin(CheatRecoShowerSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
         	  larsim_MCCheater_ParticleInventoryService_service
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandrizzleAlg    
                #dune_FDSelections_pandizzle
             )

simple_plugin(CheatRecoPfoShowerSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
         	  larsim_MCCheater_ParticleInventoryService_service
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

                   PandrizzleAlg    
                #dune_FDSelections_pandizzle
             )


simple_plugin(MuonMVARecoTrackSelector "tool"
    	  	        lardataobj_RecoBase
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas

             )

#simple_plugin(HighestEnergyRecoShowerSelector "tool"
#    	  	        lardataobj_RecoBase
#                larcorealg_Geometry
#                larcore_Geometry_Geometry_service
#                larreco_RecoAlg
#                ${ROOT_BASIC_LIB_LIST}
#                ${ART_FRAMEWORK_PRINCIPAL}
#                ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                ${ART_FRAMEWORK_CORE}
# 	  	           art_Persistency_Common canvas
# 	  	           art_Persistency_Provenance canvas
#
# 	  	           art_Utilities canvas
#             )

simple_plugin(HighestEnergyRecoVertexShowerSelector "tool"
    	  	        lardataobj_RecoBase
                larcorealg_Geometry
                larcore_Geometry_Geometry_service
                larreco_RecoAlg
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_CORE}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas

 	  	           art_Utilities canvas
                dune_AnaUtils
             )
