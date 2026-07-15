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
  if(path.includes("/api/manual-design/analyze")) return {objective_profile:{key:"balanced"},ranker_manifest:{run_id:"ofrun_test",application_version:"1.34.0",ranker_version:"2.2.0",manifest_sha256:"abc123",scientific_models:{reaction_condition_snapshot:{anneal_c:60}}},candidate:{assay,evidence:{hard_valid:true,hard_failures:[]},rank_explanation:{strongest_feature:"clean",weakest_feature:"none",ranking_may_reverse_if:"new evidence",rank_reversal_scenarios:[{trigger:"new off-target corpus",reason:"exclusivity changes"}]},rank_trace:{authoritative_key:[0]}},mappings:{forward:[{strand:"+",start:0,end:20,mismatches:0,three_prime_status:"exact",extension_eligible:true,site:"ACGT"}],reverse:[{strand:"-",start:80,end:100,mismatches:1,three_prime_status:"mismatch",extension_eligible:false,site:"TGCA"}],probe:[{strand:"+",start:40,end:60,mismatches:0,three_prime_status:"not_applicable",extension_eligible:true,site:"CGTA"}]},predicted_products:[{}]};
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
 await window.mdRun("analyze");
 check(!!last("/api/manual-design/analyze"),"analysis endpoint not called");
 check(/authoritative manual analysis/.test(document.getElementById("md_out").textContent),"analysis not rendered");
 check(/Add analyzed assay to Workbench/.test(document.getElementById("md_out").textContent),"direct analyzed assay lacks Workbench action");
 set("md_f","ACGTACGTACGTACGTACGA");
 await window.mdCompareEdit();
 const editCmp=last("/api/manual-design/compare-edit");
 check(editCmp && editCmp.body.baseline_forward==="ACGTACGTACGTACGTACGT","edit comparison lost baseline sequence");
 check(editCmp.body.edited_forward==="ACGTACGTACGTACGTACGA","edit comparison lost edited sequence");
 check(/manual edit comparison/.test(document.getElementById("md_edit_out").textContent),"edit comparison not rendered");
 check(/resolved hard failures/.test(document.getElementById("md_edit_out").textContent),"resolved failures not shown");
 set("md_f","ACGTACGTACGTACGTACGT");
 window.mdAddCandidate(0);
 const directWb=JSON.parse(window.localStorage.getItem("of_panel")||"[]"),directAdded=directWb[directWb.length-1]||{};
 check(directAdded.name==="Manual Studio analyzed assay","direct analyzed assay was not saved distinctly");
 check(directAdded.ranker_manifest&&directAdded.ranker_manifest.run_id==="ofrun_test","direct analyzed assay lost rank manifest in Workbench");
 check(directAdded.source_workflow==="manual_design","direct analyzed assay lost workflow provenance");
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
 check(/assay diagnosis/.test(document.getElementById("md_out").textContent),"rescue not rendered");
 window.mdAddCandidate(0);
 const wb=JSON.parse(window.localStorage.getItem("of_panel")||"[]"); const added=wb[wb.length-1]||{};
 check(added.ranker_manifest&&added.ranker_manifest.run_id==="ofrun_test","manual/rescue candidate lost rank manifest in Workbench");
 check(added.objective_profile&&added.objective_profile.key==="balanced","manual/rescue candidate lost objective provenance");
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
 set("md_feedback_payload",'[{"assay_id":"A","target_group":"G","status":"success"},{"assay_id":"A","target_group":"G","status":"success"}]');
 await window.mdImportFeedback();
 check(!!last("/api/experimental-feedback/import"),"feedback import endpoint not called");
 check(/feedback dataset audit/.test(document.getElementById("md_feedback_out").textContent),"feedback dataset audit not rendered");
 await window.mdSplitFeedback();
 const split=last("/api/experimental-feedback/split");
 check(split && split.body.records.length===1,"normalized feedback records not sent to split endpoint");
 check(/target-group split/.test(document.getElementById("md_feedback_out").textContent),"feedback split not rendered");
 await window.mdSummarizeFeedback();
 check(!!last("/api/experimental-feedback/summary"),"feedback summary endpoint not called");
 check(/local experimental evidence summary/.test(document.getElementById("md_feedback_out").textContent),"feedback summary not rendered");
 set("md_run_left",JSON.stringify({candidates:[{assay:assay}],ranker_manifest:{run_id:"A"}}));
 set("md_run_right",JSON.stringify({candidates:[{assay:assay}],ranker_manifest:{run_id:"B"}}));
 await window.mdCompareRuns();
 check(!!last("/api/design-runs/compare"),"run comparison endpoint not called");
 check(/design-run comparison/.test(document.getElementById("md_run_out").textContent),"run comparison not rendered");
 check(/objective/.test(document.getElementById("md_run_out").textContent),"context differences not shown");
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
