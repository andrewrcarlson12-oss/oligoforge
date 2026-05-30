const fs=require("fs"), vm=require("vm");
const html=fs.readFileSync("static/index.html","utf8");
const script=html.match(/<script>([\s\S]*)<\/script>/)[1];

// ---- minimal DOM stub good enough to run the render handlers ----
const VAL={q_ng:"10",q_len:"173",q_fold:"10",q_pts:"6",p_amp:"93",ep_max:"3000",
  au_n:"6",cons_mi:"0.6",cons_tn:"8",cons_on:"6",i_s:"225",i_e:"317",
  i_gene:"HMBS",i_org:"Aphelocoma coerulescens",i_acc:"XM_1.1",email:"x@y.com",
  b_mode:"remote",ep_mode:"remote",au_prof:"parasite_sybr",qc_role:"primer"};
const cache={};
function row(){return{querySelectorAll:()=>[{value:"HMBS_F"},{value:"GAGCTATACCCCGACCTCTG"}],
  querySelector:()=>({value:"primer"}),remove(){},parentNode:{remove(){}}};}
function el(id){if(cache[id])return cache[id];
  const e={_id:id,value:(VAL[id]!==undefined?VAL[id]:"ACGTACGTACGTACGTACGT"),checked:false,
    innerHTML:"",textContent:"",style:{},onclick:null,
    classList:{add(){},remove(){},contains:()=>false},appendChild(){},addEventListener(){},click(){},
    querySelectorAll:()=>[row(),row()],querySelector:()=>({value:"primer"})};
  cache[id]=e;return e;}
const ctx={
  document:{getElementById:el,createElement:()=>el("c"+Math.random()),
    querySelectorAll:(s)=>(s==="nav a"||s==="section")?[el("n1"),el("n2")]:[row(),row()]},
  window:{},console,Blob:function(){},URL:{createObjectURL:()=>"blob:x"},
  fetch:async()=>({ok:true,status:200,json:async()=>({})}),
  setTimeout,JSON,Math,Object,Array,Date,parseInt,parseFloat,isFinite,Number,String,Boolean
};
vm.createContext(ctx);
try{vm.runInContext(script,ctx);}catch(e){console.log("SCRIPT LOAD ERROR:",e.message);process.exit(1)}

// ---- realistic responses, incl. edge cases the endpoints really return ----
const lint=[{rule:"Tm",status:"PASS",detail:"61.8"}];
const R={
 qc:{seq:"ACGT",length:20,gc:50,tm:61.8,hairpin_dg:-0.5,hairpin_tm:30,self_dimer:-3,max_run:3,last5_gc:2,revcomp:"ACGT",note:"contains degenerate bases — approximate",lint},
 pair:{f_tm:61.8,r_tm:61.5,pair_gap:0.3,fxr:-3.2,f_self:-4,r_self:-5,lint},
 matrix:{names:["A","B"],cells:[{a:"A",b:"A",dg:-4},{a:"A",b:"B",dg:-2},{a:"B",b:"B",dg:-3}]},
 design_taq:{forward:"GAGCTATACCCCGACCTCTG",reverse:"CTTCTCTCCAATCTTGGAAAGCG",probe:"ATCTTGTCCCCAGTTGTTGACATGGCC",amplicon:93,pair_tm_gap:0.4,f_tm:61.8,r_tm:61.5,probe_tm:69,probe_offset:7.5,probe_hairpin:-1,probe_dimer_f:-3,probe_dimer_r:-4,gblock:"ACGT",f_xy:[20,40],r_xy:[100,123],profile:"IDT PrimeTime"},
 design_sybr:{forward:"GGGTTATGTATTACCTTGGACTAG",reverse:"TGGATATCTTGTAAGTGACCCA",probe:null,amplicon:122,pair_tm_gap:0.1,f_tm:55,r_tm:55,probe_tm:null,probe_offset:null,probe_hairpin:null,probe_dimer_f:null,probe_dimer_r:null,gblock:"ACGT",f_xy:[0,24],r_xy:[100,122],profile:"SYBR"},
 fetch:{query:"x",records:[{id:"XM_1.1",desc:"HMBS mRNA predicted",length:1182,seq:"ACGTACGT"}],common_region:"ACGTACGT"},
 fetch_acc:{records:[{id:"XM_1.1",desc:"HMBS",length:1182,seq:"ACGT"}]},
 intron_ok:{ok:true,verdict:"spans an exon-exon junction",info:"XM: 7 exons",junctions:[182,255,305],spanned:[255,305]},
 intron_null:{ok:null,info:"gene not found",junctions:null,spanned:[]},
 intron_false:{ok:false,verdict:"lies WITHIN a single exon",info:"info",junctions:[500],spanned:[]},
 blast_remote:{query_len:20,n_hits:1,hits:[{title:"X",identity:20,align_len:20,evalue:1e-5}]},
 blast_local:{query_len:20,n_hits:1,hits:[["acc",99,20,"1e-5","title"]]},
 blast_err:{error:"blastn not on PATH"},
 copies:{copies_per_ul:5.3e10,ng_per_ul:10,length_bp:173,series:[{point:0,copies_per_ul:5.3e10},{point:1,copies_per_ul:5.3e9}]},
 batch:{results:[{name:"HMBS",ok:true,forward:"A",reverse:"C",probe:"G",amplicon:93,f_tm:61,r_tm:61,probe_tm:69,gblock:"A"},{name:"X",ok:false,error:"no clean assay"}]},
 batch_sybr:{results:[{name:"S",ok:true,forward:"A",reverse:"C",probe:null,amplicon:122,f_tm:55,r_tm:55,probe_tm:null,gblock:"A"}]},
 order:{oligo_csv:"Name,Sequence\nHMBS_P,/56-FAM/ACGT/3IABkFQ/,100nm,HPLC",gblock_fasta:">gb\nACGT"},
 order_empty:{oligo_csv:"",gblock_fasta:""},
 cons_full:{Plas_P:{conservation:{n_placed:8,n_input:8,mean_ident:99.4,min_pct_match:90,worst_3prime:100,per_pos:[{pos:1,oligo:"C",pct_match:100,major:"C"},{pos:2,oligo:"T",pct_match:88,major:"T"}]},discrimination:{n:6,median_ident:82.5,max_ident:85,min_mismatch:3,min_3prime_mismatch:0,rows:[]}}},
 cons_nooff:{Plas_F:{conservation:{n_placed:8,n_input:8,mean_ident:93,min_pct_match:60,worst_3prime:100,per_pos:[{pos:1,oligo:"T",pct_match:100,major:"T"}]}}},
 cons_err:{error:"boom"},
 epcr_ok:{forward_hits:2,reverse_hits:2,n_products:1,products:[{subject:"chr1",size:157,left:"F",right:"R",span:[100,256]}]},
 epcr_empty:{forward_hits:1,reverse_hits:1,n_products:0,products:[]},
 epcr_err:{error:"blastn not on PATH"},
 std_ok:{slope:-3.32,intercept:40,efficiency:1,efficiency_pct:99.9,r2:1,amp_factor:2,n_points:6,dynamic_range:[10,1e6],lod_practical:100,efficiency_ok:true,r2_ok:true,slope_ok:true,levels:[{quantity:1e6,n:1,detected:1,mean_cq:20,sd_cq:null,detection_rate:100}],notes:"MIQE"},
 std_err:{error:"need >=2 detected points"},
 auto_full:{n_targets:6,n_offs:6,reference_len:478,n_candidates:1,candidates:[{score:7.5,assay:{forward:"ATGGG",reverse:"TGGAT",probe:null,amplicon:122,f_tm:55,r_tm:55,probe_info:null},conservation:{F:{mean_ident:100,worst_3prime:100},R:{mean_ident:96,worst_3prime:100}},discrimination:null}]},
 auto_probe:{n_targets:6,n_offs:0,reference_len:478,n_candidates:1,candidates:[{score:80,assay:{forward:"A",reverse:"C",probe:"G",amplicon:100,f_tm:60,r_tm:60,probe_info:{tm:69}},conservation:{F:{mean_ident:100,worst_3prime:100},R:{mean_ident:100,worst_3prime:100},P:{mean_ident:99,worst_3prime:100}},discrimination:null}]},
 auto_err:{error:"NCBI returned nothing"},
 pairspec:{forward:{n_hits:1,hits:[{title:"X",evalue:1e-5}]},reverse:{error:"blastn not on PATH"}},
};
const cases=[
 ["doQC",R.qc],["doPair",R.pair],
 ["doMatrix",R.matrix],["doMatrix(error)",R.matrix,"matrix",{error:"network drop"}],
 ["doDesign(taqman)",R.design_taq],["doDesign(sybr)",R.design_sybr],["doDesign(error)",R.design_taq,"d",{error:"no clean assay"}],
 ["doFetch",R.fetch],["doFetch(acc)",R.fetch_acc],
 ["doIntron(ok)",R.intron_ok],["doIntron(null)",R.intron_null],["doIntron(false)",R.intron_false],
 ["doBlast(remote)",R.blast_remote],["doBlast(local)",R.blast_local],["doBlast(err)",R.blast_err],
 ["doCopies",R.copies],
 ["doBatch",R.batch],["doBatch(sybr)",R.batch_sybr],
 ["doOrder",R.order],["doOrder(empty)",R.order_empty],
 ["doCons(full)",R.cons_full],["doCons(nooff)",R.cons_nooff],["doCons(error)",R.cons_err],
 ["doEpcr(ok)",R.epcr_ok],["doEpcr(empty)",R.epcr_empty],["doEpcr(err)",R.epcr_err],
 ["doStd(ok)",R.std_ok],["doStd(err)",R.std_err],
 ["doAuto(full)",R.auto_full],["doAuto(probe)",R.auto_probe],["doAuto(error)",R.auto_err],
 ["doPairSpec",R.pairspec],
];
const fnOf={doMatrix:"doMatrix",doDesign:"doDesign",doFetch:"doFetch",doIntron:"doIntron",doBlast:"doBlast",doCons:"doCons",doEpcr:"doEpcr",doStd:"doStd",doAuto:"doAuto",doQC:"doQC",doPair:"doPair",doCopies:"doCopies",doBatch:"doBatch",doOrder:"doOrder"};
(async()=>{
 let pass=0,fail=[];
 for(const [label,resp] of cases){
   const fn=label.replace(/\(.*/,"");
   ctx.api=async()=>resp;
   try{ const r=ctx[fn]==null?(()=>{throw new Error("handler "+fn+" not defined")})():await ctx[fn](fn==="doPairSpec"?"AAA":undefined, "CCC"); pass++; }
   catch(e){ fail.push(label+" -> "+e.message); }
 }
 console.log(pass+" handler renders OK, "+fail.length+" threw");
 fail.forEach(f=>console.log("  ✗ "+f));
 if(!fail.length)console.log("  (no runtime render crashes)");
})();
