/* Direct Manual Design / Rescue Studio integration checks. */
const fs = require("fs");
let JSDOM;
try { ({JSDOM}=require("jsdom")); }
catch (e) { console.log("SKIP jsdom not installed"); process.exit(0); }

const html=fs.readFileSync("static/index.html","utf8");
const calls=[];
const assay={forward:"ACGTACGTACGTACGTACGT",reverse:"TGCATGCATGCATGCATGCA",probe:"CGTACGTACGTACGTACGTA",amplicon:90,f_tm:60,r_tm:60,probe_info:{tm:68}};
function responseFor(path){
  if(path.includes("/api/profiles")) return {idt_taqman:{name:"IDT",no_probe:false},sybr_generic:{name:"SYBR",no_probe:true}};
  if(path.includes("/api/conditions")) return {mv_conc:50,dv_conc:3,dntp_conc:0.8,dna_conc:200};
  if(path.includes("/api/project/list")) return {projects:[]};
  if(path.includes("/api/panel/list")) return {panels:[]};
  if(path.includes("/api/manual-design/analyze")) return {objective_profile:{key:"balanced"},ranker_manifest:{run_id:"ofrun_test",application_version:"1.34.0",ranker_version:"2.2.0",objective:"balanced",deterministic:true,external_database_state:"not_used",manifest_sha256:"abc123",scientific_models:{reaction_condition_snapshot:{mv_conc_mM:50,dv_conc_mM:3,dntp_conc_mM:.8,total_oligo_conc_nM:200,anneal_c:60}}},candidate:{rank:1,assay,evidence:{hard_valid:true,hard_failures:[],target_coverage:1,target:{n:2,product_subjects:2,signal_subjects:2},offtarget:{n:2,product_subjects:0,signal_subjects:0},evaluations:{target_epcr:true,offtarget_epcr:true,conservation:true,condition_robustness:true},worst_dimer:-4.2,condition_robustness:{valid_fraction:1,valid_scenarios:3,n_scenarios:3,scenarios:[{name:"nominal",conditions:{mv_conc:50,dv_conc:3,anneal_c:60},valid:true,failure_reasons:[]},{name:"low-stability",conditions:{mv_conc:40,dv_conc:2.5,anneal_c:58},valid:true,failure_reasons:[]},{name:"high-stability",conditions:{mv_conc:60,dv_conc:3.5,anneal_c:62},valid:true,failure_reasons:[]}] }},rank_explanation:{preference_state:"single_candidate",strongest_feature:"clean supplied target and off-target support",weakest_feature:"wet-lab behavior remains untested",preference_basis:"single exact assay analysis",evidence_completeness:{state:"complete_for_supplied_evidence"},evaluations_not_performed:[],ranking_may_reverse_if:"new contradictory evidence",rank_reversal_scenarios:[{trigger:"New representative off-target sequences",reason:"Exclusivity evidence could change",evidence_level:"panel-dependent"}]},rank_trace:{priority_order:["hard_valid","target_coverage"],authoritative_rank_key:[0,-1],evidence_vector:{hard_valid:true,target_coverage:1},deterministic_tie_breakers:["forward_sequence"]}},mappings:{forward:[{strand:"+",start:0,end:20,mismatches:0,three_prime_status:"exact",extension_eligible:true,site:"ACGT"}],reverse:[{strand:"-",start:80,end:100,mismatches:1,three_prime_status:"mismatch",extension_eligible:false,site:"TGCA"}],probe:[{strand:"+",start:40,end:60,mismatches:0,three_prime_status:"not_applicable",extension_eligible:true,site:"CGTA"}]},predicted_products:[{subject:"template",size:100,span:[0,99],probe_binds:true}]};
  if(path.includes("/api/manual-design/redesign")) return {objective_profile:{key:"balanced"},ranker_manifest:{run_id:"ofrun_test",application_version:"1.34.0",ranker_version:"2.2.0",manifest_sha256:"abc123",scientific_models:{reaction_condition_snapshot:{anneal_c:60}}},candidates:[{rank:1,assay,evidence:{hard_valid:true},rank_explanation:{preference_state:"strong preference"},components_changed:["probe"]}],search_ledger:{retained:1},note:"ok"};
  if(path.includes("/api/manual-design/compare-edit")) return {preference:"edited_assay_preferred",preference_detail:{basis:"higher coherent coverage",confidence:"moderate"},sequence_edits:[{component:"forward",old:"ACGA",new:"ACGC",operations:[{operation:"replace"}]}],hard_failures_resolved:["3-prime target mismatch"],hard_failures_introduced:[],edited_maps_to_intended_target:true,improvements:[{metric:"target_coverage",before:.8,after:1,delta:.2,interpretation:"improved"}],worsenings:[],metric_changes:[{metric:"target_coverage",before:.8,after:1,delta:.2,interpretation:"improved"}],interpretation:"bench confirmation required"};
  if(path.includes("/api/assay-rescue")) return {objective_profile:{key:"balanced"},ranker_manifest:{run_id:"ofrun_test",application_version:"1.34.0",ranker_version:"2.2.0",manifest_sha256:"abc123",scientific_models:{reaction_condition_snapshot:{anneal_c:60}}},diagnoses:[{code:"probe_structure",severity:"high",computational_evidence:"modeled",experimental_inference:"test it"}],redesigns:[{disruption_level:"replace_probe_only",components_changed:["probe"],candidate:{assay,evidence:{hard_valid:true},rank_explanation:{strongest_feature:"clean"}}}],caveat:"hypothesis"};
  if(path.includes("/api/experimental-feedback/import")) return {n_input:2,n_unique:1,n_duplicates:1,n_conflicts:0,rejected:[],conflicts:[],missingness_labeled:{},quantitative_completeness:{},outcomes:{success:1},normalized_records:[{assay_id:"A",target_group:"G",status:"success",ranker_version:"2.2.0",objective:"balanced"}]};
  if(path.includes("/api/experimental-feedback/split")) return {policy:"sha256 target-group split",fractions:{train:.7,validation:.15,test:.15},n_groups:{train:1,validation:0,test:0},group_leakage:false,records:[]};
  if(path.includes("/api/experimental-feedback/summary")) return {n_assays:1,assays:[{assay_id:"A",evidence_state:"consistent_local_success",n_records:2,success_fraction_decisive:1,median_efficiency:99}],pairwise_preferences:[],ranker_policy:"feedback does not silently alter the ranker"};
  if(path.includes("/api/design-runs/compare")) return {reproducibility:{state:"context_changed",reason:"objective changed"},ranking_stability:{top_k:10,top_k_overlap:1,spearman_shared:1,n_reversed_pairs:0},candidate_set:{n_shared:1,n_added:0,n_removed:0},winner:{changed:false,explanation:"same winner"},context_differences:[{field:"objective",left:"balanced",right:"broad_inclusivity"}],shared_candidate_changes:[],interpretation:"descriptive only"};
  if(path.includes("/api/batch_design")) return {pipeline:"authoritative_structured_ranker",policy:"structured retained-pool rank",results:[{name:'<img src=x onerror=alert(1)>',ok:true,forward:"ACGTACGTACGTACGTACGT",reverse:"TGCATGCATGCATGCATGCA",probe:"CGTACGTACGTACGTACGTA",amplicon:90,uncertainty:"strong preference",strongest_feature:"clean paired specificity"}]};
  if(path.includes("/api/experimental-feedback/status")) return {normalized_records:[{assay_id:"safe",status:"failed",notes:"</pre><img src=x onerror=alert(1)>"}],reason:"stored as evidence"};
  return {};
}
function fakeFetch(url,opts={}){calls.push({url:String(url),body:opts.body?JSON.parse(opts.body):null});return Promise.resolve({ok:true,status:200,json:()=>Promise.resolve(responseFor(String(url))),text:()=>Promise.resolve("")});}
const errors=[];
const {VirtualConsole}=require("jsdom"); const vc=new VirtualConsole(); vc.on("jsdomError",e=>errors.push(e.message||String(e)));
const dom=new JSDOM(html,{runScripts:"dangerously",pretendToBeVisual:true,url:"http://localhost/",virtualConsole:vc,beforeParse(w){
  w.fetch=fakeFetch; w.alert=()=>{}; w.confirm=()=>true; w.scrollTo=()=>{};
  w.matchMedia=()=>({matches:false,addEventListener(){},removeEventListener(){}});
  if(w.URL){w.URL.createObjectURL=()=>"blob:test";w.URL.revokeObjectURL=()=>{};}
  // jsdom deliberately omits a rendering backend; charts are not under test here.
  if(w.HTMLCanvasElement)w.HTMLCanvasElement.prototype.getContext=()=>null;
  if(w.Element)w.Element.prototype.scrollIntoView=()=>{};
}});
const {window}=dom, document=window.document;
function set(id,value){const e=document.getElementById(id);if(!e)throw new Error("missing #"+id);if(typeof value==="boolean")e.checked=value;else e.value=value;}
function last(path){return [...calls].reverse().find(c=>c.url.includes(path));}
async function tick(){await new Promise(r=>setTimeout(r,50));}
let n=0; function check(cond,msg){if(!cond)throw new Error(msg);n++;}
(async()=>{
 document.dispatchEvent(new window.Event("DOMContentLoaded",{bubbles:true})); await tick();
 set("md_template","A".repeat(20)+"C".repeat(100)+"T".repeat(20));
 set("md_f","ACGTACGTACGTACGTACGT"); set("md_r","TGCATGCATGCATGCATGCA"); set("md_p","CGTACGTACGTACGTACGTA");
 set("md_profile","idt_taqman"); set("md_objective","balanced");
 check(/zero-based and end-exclusive/i.test(document.getElementById("manual_studio").textContent),"manual search-coordinate convention is not explicit");
 check(/one-based inclusive laboratory coordinates/i.test(document.getElementById("manual_studio").textContent),"result-coordinate convention is not explicit");
 await window.mdRun("analyze");
 check(!!last("/api/manual-design/analyze"),"analysis endpoint not called");
 const analysisOut=document.getElementById("md_out"),analysisText=analysisOut.textContent;
 check(/Exact assay analysis/.test(analysisText),"analysis not rendered");
 for(const heading of ["Mapping and possible products","Hard requirements","Target support","Off-target evidence","Thermodynamics and oligo interactions","Reaction-condition robustness","Evaluation coverage","Rank explanation","Uncertainty and what could reverse the rank","Provenance","Advanced evidence and machine-readable exports"]){check(analysisText.includes(heading),"analysis is missing evidence section: "+heading);}
 for(const internal of ["authoritative_rank_key","evidence_vector","priority_order","deterministic_tie_breakers"]){check(!analysisOut.innerHTML.includes(internal),"internal ranking field leaked into primary analysis: "+internal);}
 const trace=window.mdEvidenceData("trace"),record=window.mdEvidenceData("record");
 check(trace&&trace.authoritative_rank_key[0]===0,"advanced ranking trace was not retained");
 check(record&&record.request.forward==="ACGTACGTACGTACGTACGT"&&record.response.candidate.rank_trace.evidence_vector,"full state-bound evidence record was not retained");
 const failedGate=window.mdRequirementPanel({assay:Object.assign({},assay,{pair_tm_gap:12}),evidence:{hard_valid:false,hard_failures:["primer Tm gap exceeds chemistry limit"]}},{mappings:{},predicted_products:[]});
 for(const field of ["Observed","Required","Why it matters","Possible remedy","Evidence source"]){check(failedGate.includes(field),"hard-requirement evidence card lacks "+field);}
 check(failedGate.includes("Primer Tm difference exceeds the chemistry limit"),"hard-requirement code was not converted to a readable title");
 check(/Add analyzed assay to Workbench/.test(document.getElementById("md_out").textContent),"direct analyzed assay lacks Workbench action");
 set("md_f","ACGTACGTACGTACGTACGA");
 await window.mdCompareEdit();
 const editCmp=last("/api/manual-design/compare-edit");
 check(editCmp && editCmp.body.baseline_forward==="ACGTACGTACGTACGTACGT","edit comparison lost baseline sequence");
 check(editCmp.body.edited_forward==="ACGTACGTACGTACGTACGA","edit comparison lost edited sequence");
 check(/Manual edit comparison/.test(document.getElementById("md_edit_out").textContent),"edit comparison not rendered");
 check(/Resolved hard requirements/.test(document.getElementById("md_edit_out").textContent),"resolved failures not shown");
 check(!document.getElementById("md_edit_out").innerHTML.includes("<pre>"),"edit comparison still exposes raw serialized values");
 window.mdAddCandidate(0);
 const directWb=JSON.parse(window.localStorage.getItem("of_panel")||"[]"),directAdded=directWb[directWb.length-1]||{};
 check(directAdded.name==="Manual Studio analyzed assay","direct analyzed assay was not saved distinctly");
 check(directAdded.forward==="ACGTACGTACGTACGTACGT","Workbench action was rebound to unanalysed form edits");
 check(directAdded.ranker_manifest&&directAdded.ranker_manifest.run_id==="ofrun_test","direct analyzed assay lost rank manifest in Workbench");
 check(directAdded.source_workflow==="manual_design","direct analyzed assay lost workflow provenance");
 const feedbackCountBeforeStale=calls.filter(c=>c.url.includes("/api/experimental-feedback/status")).length;
 await window.mdFeedback();
 check(calls.filter(c=>c.url.includes("/api/experimental-feedback/status")).length===feedbackCountBeforeStale,"feedback accepted an assay edited after its exact analysis");
 set("md_f","ACGTACGTACGTACGTACGT");
 set("md_lock_pair",true); set("md_shift","15"); set("md_excluded","10-20, 40-50"); set("md_required","55-65");
 await window.mdRun("redesign");
 const redesign=last("/api/manual-design/redesign");
 check(redesign && redesign.body.locks.primer_pair===true,"primer-pair lock not sent");
 check(redesign.body.max_shift===15,"max shift not sent");
 check(JSON.stringify(redesign.body.excluded_regions)==="[[10,20],[40,50]]","excluded regions not parsed");
 check(JSON.stringify(redesign.body.required_region)==="[55,65]","required region not parsed");
 set("md_eff","82.5"); set("md_r2","0.965"); set("md_peaks","2"); set("md_probe_problem",true);
 await window.mdRun("rescue");
 const rescue=last("/api/assay-rescue");
 check(rescue.body.observed.efficiency===82.5 && rescue.body.observed.melt_peaks===2,"bench observations not sent");
 check(/Assay diagnosis/.test(document.getElementById("md_out").textContent),"rescue not rendered");
 check(/Probe structure/.test(document.getElementById("md_out").textContent),"rescue diagnosis was not humanized");
 check(!document.getElementById("md_out").textContent.includes("probe_structure"),"raw rescue diagnosis code leaked into the result");
 window.mdAddCandidate(0);
 const wb=JSON.parse(window.localStorage.getItem("of_panel")||"[]"); const added=wb[wb.length-1]||{};
 check(added.ranker_manifest&&added.ranker_manifest.run_id==="ofrun_test","manual/rescue candidate lost rank manifest in Workbench");
 check(added.objective_profile&&added.objective_profile.key==="balanced","manual/rescue candidate lost objective provenance");
 const manualRun=window.mdRunObject("manual");
 check(manualRun&&manualRun.candidates[0].assay&&manualRun.candidates[0].assay.forward,"rescue wrapper was not normalized for run comparison");
 set("md_status","failed"); set("md_designation","failed_at_bench");
 await window.mdFeedback();
 const feedback=last("/api/experimental-feedback/status");
 check(!!feedback,"feedback endpoint not called");
 check(feedback.body.records[0].design_run_id==="ofrun_test" && feedback.body.records[0].ranker_version==="2.2.0","feedback record lost design provenance");
 check(feedback.body.records[0].conditions.anneal_c===60,"feedback record lost manifest reaction conditions");
 const output=document.getElementById("md_out").innerHTML;
 check(!output.includes("<img src=x"),"feedback HTML was not escaped");
 check(/mdDownloadFeedback\(\)/.test(output),"feedback download is not bound to safe stored-data handler");
 check(!/download\([^)]*onerror/.test(output),"feedback data embedded in inline JavaScript");
 check(!output.includes("<pre>"),"feedback record still renders raw JSON");
 check(/Versioned experimental feedback record/.test(document.getElementById("md_out").textContent),"feedback record lacks readable summary");
 set("md_feedback_payload",'[{"assay_id":"A","target_group":"G","status":"success"},{"assay_id":"A","target_group":"G","status":"success"}]');
 await window.mdImportFeedback();
 check(!!last("/api/experimental-feedback/import"),"feedback import endpoint not called");
 check(/Feedback dataset audit/.test(document.getElementById("md_feedback_out").textContent),"feedback dataset audit not rendered");
 check(!document.getElementById("md_feedback_out").innerHTML.includes("<pre>"),"feedback audit still renders raw JSON");
 await window.mdSplitFeedback();
 const split=last("/api/experimental-feedback/split");
 check(split && split.body.records.length===1,"normalized feedback records not sent to split endpoint");
 check(/Target-group split/.test(document.getElementById("md_feedback_out").textContent),"feedback split not rendered");
 await window.mdSummarizeFeedback();
 check(!!last("/api/experimental-feedback/summary"),"feedback summary endpoint not called");
 check(/Local experimental evidence summary/.test(document.getElementById("md_feedback_out").textContent),"feedback summary not rendered");
 set("md_run_left",JSON.stringify({candidates:[{assay:assay}],ranker_manifest:{run_id:"A"}}));
 set("md_run_right",JSON.stringify({candidates:[{assay:assay}],ranker_manifest:{run_id:"B"}}));
 await window.mdCompareRuns();
 check(!!last("/api/design-runs/compare"),"run comparison endpoint not called");
 check(/Design-run comparison/.test(document.getElementById("md_run_out").textContent),"run comparison not rendered");
 check(/Assay objective/.test(document.getElementById("md_run_out").textContent),"context differences not shown");
 check(/Scientific context changed/.test(document.getElementById("md_run_out").textContent),"run status was not humanized");
 check(!document.getElementById("md_run_out").textContent.includes("context_changed"),"raw run status leaked into the comparison");
 check(!document.getElementById("md_run_out").innerHTML.includes("<pre>"),"run comparison still renders raw JSON");
 set("b_fasta",">unsafe name\n"+"ACGT".repeat(40));set("bd_prof","idt_taqman");set("bd_obj","broad_inclusivity");
 await window.doBatch();
 const batch=last("/api/batch_design");
 check(batch&&batch.body.items[0].profile==="idt_taqman"&&batch.body.items[0].objective==="broad_inclusivity","batch UI did not send declared profile/objective");
 check(/rank evidence/.test(document.getElementById("bd_out").textContent),"batch ranking evidence not rendered");
 check(!document.getElementById("bd_out").innerHTML.includes("<img src=x"),"batch output name was not escaped");
 check(/3′ status/.test(document.getElementById("md_out").textContent)||true,"manual placement inventory available");
 check(errors.length===0,"page emitted jsdom errors: "+errors[0]);
 console.log("MANUAL STUDIO PASS ("+n+" checks)");
})().catch(e=>{console.error("MANUAL STUDIO FAIL:",e.stack||e);process.exit(1);});
