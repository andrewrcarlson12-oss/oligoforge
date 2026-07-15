const fs=require("fs"),vm=require("vm");
const script=fs.readFileSync("static/index.html","utf8").match(/<script>([\s\S]*)<\/script>/)[1];
const VAL={pj_name:"Smith frog panel",pj_sel:"Smith frog panel",cond_dv:"3",cond_mv:"50",cond_dntp:"0.8",cond_dna:"200",email:"x@y.com"};
const store={}; const cache={};
function el(id){if(cache[id])return cache[id];cache[id]={_id:id,value:(VAL[id]!==undefined?VAL[id]:"ACGT"),checked:false,innerHTML:"",textContent:"",className:"",style:{},classList:{add(){},remove(){},contains:()=>false},appendChild(){},addEventListener(){},click(){},querySelectorAll:()=>[],querySelector:()=>({value:"primer"})};return cache[id];}
const ctx={document:{getElementById:el,createElement:()=>el("c"+Math.random()),querySelectorAll:(s)=>(s==="nav a"||s==="section")?Array.from({length:17},(_,i)=>el("n"+i)):[]},window:{},console,Blob:function(){},URL:{createObjectURL:()=>"x"},fetch:async()=>({ok:true,status:200,json:async()=>({})}),localStorage:{getItem:k=>store[k]||null,setItem:(k,v)=>{store[k]=v;}},alert:()=>{},confirm:()=>true,setTimeout:(f)=>0,JSON,Math,Object,Array,Date,parseInt,parseFloat,isFinite,Number,String,Boolean};
vm.createContext(ctx); vm.runInContext(script,ctx);
const RESP={"/api/conditions":{mv_conc:50,dv_conc:3,dntp_conc:0.8,dna_conc:200},
 "/api/project/list":{projects:[{name:"Smith frog panel",n:5,saved:"2026-05-30 16:00"}]},
 "/api/project/save":{saved:"Smith frog panel",n:5},
 "/api/project/load":{name:"Smith frog panel",assays:[{name:"IL1B",gene:"IL1B",organism:"Lithobates pipiens",forward:"ACGT",reverse:"TGCA",probe:null,amplicon:110,chemistry:"SYBR"}]},
 "/api/project/delete":{deleted:"Smith frog panel"}};
ctx.api=async(p)=>RESP[p]||{};
(async()=>{let ok=0,bad=[];const T=async(l,fn)=>{try{await fn();ok++;console.log("  \u2713 "+l);}catch(e){bad.push(l);console.log("  \u2717 "+l+" -> "+e.message);}};
 await T("loadProjectList (populate select)",async()=>{await ctx.loadProjectList();if(el("pj_sel").innerHTML.indexOf("Smith frog")<0)throw new Error("list not rendered");});
 await T("seedFSJ then saveProject",async()=>{ctx.seedFSJ();await ctx.saveProject();});
 await T("loadProject (replaces workbench w/ frog panel)",async()=>{await ctx.loadProject();if(el("wk_out").innerHTML.indexOf("IL1B")<0)throw new Error("loaded panel not rendered");});
 await T("deleteProject",async()=>{await ctx.deleteProject();});
 console.log("\n"+ok+" ok, "+bad.length+" failed"+(bad.length?": "+bad.join(", "):""));
})();
