/* Real-DOM integration coverage for the visible Validation Studio and Assurance workspaces. */
const fs=require("fs");
let JSDOM;
try{({JSDOM}=require("jsdom"));}catch(e){console.log("SKIP jsdom not installed");process.exit(0);}
const html=fs.readFileSync("static/index.html","utf8"),calls=[];
function parsed(opts){try{return opts&&opts.body?JSON.parse(opts.body):null;}catch(e){return null;}}
function response(path,body){
 if(path==="/api/profiles")return {idt_taqman:{name:"IDT PrimeTime",no_probe:false},sybr_generic:{name:"SYBR",no_probe:true}};
 if(path==="/api/conditions")return {mv_conc:50,dv_conc:3,dntp_conc:0.8,dna_conc:200};
 if(path==="/api/project/list"||path==="/api/panel/list")return path.includes("project")?{projects:[]}:{panels:[]};
 if(path==="/api/autodesign/limits")return {queue_capacity:8,primary_timeout_seconds:240,blast_timeout_seconds:360,max_records_to_fetch:30,terminal_ttl_seconds:1800};
 if(path==="/healthz")return {ok:true,version:"1.37.0"};
 if(path==="/api/validation-studio/plan"){
  const c=body.candidates;
  const states={};states[c[0].candidate_id]={state:"signal_product",basis:"modeled product"};states[c[1].candidate_id]={state:"no_product",basis:"terminal-primer concern"};
  return {plan:{schema_version:"oligoforge-validation-plan/v1",plan_sha256:"planhash",objective:body.objective,candidates:c,selected_cases:[{case_id:"target-A",name:"target-A",role:"target",group:"lineage-A",source_type:"synthetic",selection_rationale:"tests candidate-specific modeled target support",distinguishes_candidate_pairs:[[c[0].candidate_id,c[1].candidate_id]],modeled_states:states}],selection_status:{method:"deterministic greedy pair-coverage with group diversity",globally_optimal:false,n_evaluated:2,n_informative:2,n_selected:1},plate_layout:{plate_format:96,n_wells:6,randomization_seed:body.seed,warnings:["2 test wells are on a plate edge"],wells:[{well:"A1",well_type:"test",candidate_id:c[0].candidate_id,case_id:"target-A",replicate:1,edge_well:true},{well:"B1",well_type:"test",candidate_id:c[1].candidate_id,case_id:"target-A",replicate:1,edge_well:false},{well:"C1",well_type:"test",candidate_id:c[0].candidate_id,case_id:"target-A",replicate:2,edge_well:false},{well:"D1",well_type:"test",candidate_id:c[1].candidate_id,case_id:"target-A",replicate:2,edge_well:false},{well:"E1",well_type:"no_template_control",candidate_id:c[0].candidate_id,case_id:"NTC",replicate:1,edge_well:false},{well:"F1",well_type:"no_template_control",candidate_id:c[1].candidate_id,case_id:"NTC",replicate:1,edge_well:false}]},limitations:["Declared cases only."]},plate_csv:"well,well_type,candidate_id,case_id,observed\nA1,test,candidate-a,target-A,\n"};
 }
 if(path==="/api/validation-studio/interpret")return {schema_version:"oligoforge-validation-interpretation/v1",controls_valid:true,invalid_control_reasons:[],conclusion_strength:"moderate for this declared experiment",conclusion:"Moderate for this declared experiment favors candidate-a; this is not a global assay-performance claim",candidate_summary:[{candidate_id:"candidate-a",supported:1,contradicted:0,missing:0},{candidate_id:"candidate-b",supported:0,contradicted:1,missing:0}],predictions_supported:[{candidate_id:"candidate-a"}],predictions_contradicted:[{candidate_id:"candidate-b",case_id:"target-A",prediction:"no_product",n_amplified:1,n_not_amplified:0,mean_cq:25}],uncertainties_remaining:[],limitations:["Interpretation applies only to this declared experiment."]};
 if(path==="/api/assurance/assaysbom")return {schema_version:"oligoforge-assaysbom/v1",assaysbom_id:"ofsbom_fixture",portfolio_version:"1",review_state:"unreviewed",content_sha256:"abcdef",scope_statement:"Molecular and computational bill of materials; not proof of assay performance.",assays:[{assay_id:"assay-A",name:"Synthetic probe assay",assay_type:"probe",components:[{role:"forward_primer",order_sequence:"ACGTACGTACGTACGTACGT",locked_legacy_component:false},{role:"reverse_primer",order_sequence:"TTTTGGGGCCCCAAAATTTT",locked_legacy_component:false},{role:"probe",order_sequence:"GATTACAGATTACAGATTAC",locked_legacy_component:false}]}]};
 if(path==="/api/assurance/snapshots"){
  const ids={"Baseline target snapshot":"ofsnap_bt","Follow-up target snapshot":"ofsnap_ct","Baseline off-target snapshot":"ofsnap_bo","Follow-up off-target snapshot":"ofsnap_co"};
  return {schema_version:"oligoforge-sequence-snapshot/v1",snapshot_id:ids[body.name],baseline_snapshot_id:body.baseline_snapshot_id||null,content_sha256:"snapsha",metrics:{raw_record_count:body.name.includes("Follow-up")?2:1,unique_sequence_count:body.name.includes("Follow-up")?2:1,exact_duplicate_count:0,rejected_record_count:0}};
 }
 if(path==="/api/assurance/snapshots/delta")return {schema_version:"oligoforge-snapshot-delta/v1",delta_sha256:"deltasha",counts:{added:1,removed:0,unchanged:1}};
 if(path==="/api/assurance/drift-scan")return {schema_version:"oligoforge-drift-scan/v1",scan_id:"ofscan_fixture",state:"Possible target dropout",action_review_recommended:true,model_version:"specificity-test",search_completeness:{status:"complete"},limitations:["Sequence-level modeled evidence only."],reason_records:[{code:"new_target_lost_coherent_product",assay_id:"assay-A",predicted_consequence:"possible target dropout under the declared product model",record:{record_id:"target-2",group:"lineage-2",affected_components:["assay-A-forward-primer-1"]}}],assay_results:[{assay_id:"assay-A",baseline_target:{signal_products:1,n_records:1},current_target:{signal_products:1,n_records:2},baseline_off_target:{signal_products:0,n_records:1},current_off_target:{signal_products:1,n_records:2}}]};
 if(path==="/api/assurance/ofvr")return {records:[{schema_version:"oligoforge-ofvr/v1",ofvr_id:"OFVR-2026-FIXTURE",reason_code:"new_target_lost_coherent_product",predicted_consequence:"possible target dropout",review_status:"unreviewed",evidence_strength:"computationally_supported_within_declared_model"}]};
 if(path==="/api/assurance/package")return {package:{schema_version:"oligoforge-assurance-evidence-package/v1",package_id:"ofpkg_fixture",package_sha256:"pkgsha",scope_statement:"Reproducible computational evidence package; expert and laboratory confirmation required.",manifest:[{category:"assaysbom",index:0,schema_version:"oligoforge-assaysbom/v1",sha256:"a"},{category:"validation_plans",index:0,schema_version:"oligoforge-validation-plan/v1",sha256:"b"}]},verification:{valid:true,artifact_checks:[{valid:true},{valid:true}]},html:"<!doctype html><h1>verified</h1>"};
 return {};
}
function fakeFetch(url,opts){
 const path=String(url),body=parsed(opts);calls.push({path:path,body:body,method:(opts&&opts.method)||"GET"});
 return Promise.resolve({ok:true,status:200,json:()=>Promise.resolve(response(path,body)),text:()=>Promise.resolve("")});
}
const errors=[];const {VirtualConsole}=require("jsdom"),vc=new VirtualConsole();vc.on("jsdomError",e=>errors.push(e.message||String(e)));
const dom=new JSDOM(html,{runScripts:"dangerously",pretendToBeVisual:true,url:"http://localhost/",virtualConsole:vc,beforeParse(w){
 w.fetch=fakeFetch;w.alert=()=>{};w.confirm=()=>true;w.scrollTo=()=>{};w.matchMedia=()=>({matches:false,addEventListener(){},removeEventListener(){}});
 if(w.URL){w.URL.createObjectURL=()=>"blob:test";w.URL.revokeObjectURL=()=>{};}if(w.HTMLCanvasElement)w.HTMLCanvasElement.prototype.getContext=()=>null;if(w.Element)w.Element.prototype.scrollIntoView=()=>{};
}});
const {window}=dom,document=window.document;
function last(path){return [...calls].reverse().find(x=>x.path===path);}
function check(value,label){if(!value)throw new Error(label);console.log("  ✓ "+label);}
async function tick(){await new Promise(r=>setTimeout(r,40));}
(async()=>{
 document.dispatchEvent(new window.Event("DOMContentLoaded",{bubbles:true}));await tick();
 check(document.getElementById("validation")&&document.getElementById("assurance"),"lifecycle pages are present");
 check([...document.querySelectorAll("#nav a")].some(x=>/Validation Studio/.test(x.textContent))&&[...document.querySelectorAll("#nav a")].some(x=>/Assurance/.test(x.textContent)),"lifecycle pages are visible in navigation");
 const lifecycleControls=[...document.querySelectorAll("#validation input, #validation textarea, #validation select, #assurance input, #assurance textarea, #assurance select")].filter(x=>x.type!=="hidden"&&!/display\s*:\s*none/i.test(x.getAttribute("style")||""));
 const unnamed=lifecycleControls.filter(x=>!(x.getAttribute("aria-label")||x.getAttribute("aria-labelledby")||(x.id&&document.querySelector('label[for="'+x.id+'"]'))||x.closest("label")));
 check(lifecycleControls.length>20&&unnamed.length===0,"every visible lifecycle form control has an accessible name");
 check(document.querySelectorAll("#validation [role=status][aria-live=polite], #assurance [role=status][aria-live=polite]").length>=9,"workflow status changes are announced politely");
 check(document.querySelectorAll("#vs_steps [aria-current=step], #as_steps [aria-current=step]").length===2,"each workflow exposes one current step");
 window.vsLoadExample();await window.vsBuildPlan();
 check(last("/api/validation-studio/plan").body.candidates.length===2,"Validation Studio submits parsed candidates");
 check(/Why these cases were selected/.test(document.getElementById("vs_out").textContent),"plan explains informative-case selection");
 check(document.querySelectorAll("#vs_out .lc-well").length>=96,"96-well plate map is rendered");
 check(!document.getElementById("vs_out").innerHTML.includes("<pre>"),"plan primary view does not dump raw JSON");
 document.getElementById("vs_results_csv").value="well,observed\nA1,amplified";await window.vsInterpret();
 check(/favors candidate-a/i.test(document.getElementById("vs_interpret_out").textContent),"completed experiment receives conservative candidate interpretation");
 window.asLoadExample();await window.asRegister(false);
 check(last("/api/assurance/assaysbom").body.assay.assays[0].forward.length===20,"displayed assay is registered as an AssaySBOM");
 await window.asCreateSnapshot("baseline_target");await window.asCreateSnapshot("current_target");
 check(last("/api/assurance/snapshots").body.baseline_snapshot_id==="ofsnap_bt","follow-up target preserves baseline linkage");
 await window.asCreateSnapshot("baseline_offtarget");await window.asCreateSnapshot("current_offtarget");
 check(last("/api/assurance/snapshots").body.baseline_snapshot_id==="ofsnap_bo","follow-up off-target preserves baseline linkage");
 check(/Target snapshot delta/.test(document.getElementById("as_snapshot_out").textContent),"snapshot deltas are visible");
 await window.asRunDrift();
 check(/Possible target dropout/.test(document.getElementById("as_drift_out").textContent),"DriftGuard state is visible");
 check(/OFVR-2026-FIXTURE/.test(document.getElementById("as_drift_out").textContent),"Molecular Vulnerability Record is visible");
 check(!document.getElementById("as_drift_out").innerHTML.includes("<pre>"),"Assurance primary view does not dump raw JSON");
 await window.asBuildPackage();
 check(last("/api/assurance/package").body.validation_plans.length===1,"evidence package includes the Validation Studio plan");
 check(/Package verified/.test(document.getElementById("as_package_out").textContent),"evidence package verification is visible");
 check(errors.length===0,"page emits no jsdom errors");
 console.log("LIFECYCLE UI PASS (20 checks)");
})().catch(e=>{console.error("LIFECYCLE UI FAIL:",e.stack||e);process.exit(1);});
