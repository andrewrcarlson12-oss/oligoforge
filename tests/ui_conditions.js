const fs=require("fs"),vm=require("vm");
const script=fs.readFileSync("static/index.html","utf8").match(/<script>([\s\S]*)<\/script>/)[1];
const VAL={std_in:"1e6, 20.1\n1e5, 23.4\n1e4, 26.7\n1e3, 30.0",email:"x@y.com",cond_dv:"3",cond_mv:"50",cond_dntp:"0.8",cond_dna:"200"};
const store={}; const cache={};
function el(id){if(cache[id])return cache[id];cache[id]={_id:id,value:(VAL[id]!==undefined?VAL[id]:"ACGTACGTACGTACGT"),checked:false,innerHTML:"",textContent:"",style:{},classList:{add(){},remove(){},contains:()=>false},appendChild(){},addEventListener(){},click(){},querySelectorAll:()=>[],querySelector:()=>({value:"primer"})};return cache[id];}
const ctx={document:{getElementById:el,createElement:()=>el("c"+Math.random()),querySelectorAll:(s)=>(s==="nav a"||s==="section")?Array.from({length:17},(_,i)=>el("n"+i)):[]},window:{},console,Blob:function(){},URL:{createObjectURL:()=>"x"},fetch:async()=>({ok:true,status:200,json:async()=>({})}),localStorage:{getItem:k=>store[k]||null,setItem:(k,v)=>{store[k]=v;}},alert:()=>{},confirm:()=>true,setTimeout:(f)=>f&&0,JSON,Math,Object,Array,Date,parseInt,parseFloat,isFinite,Number,String,Boolean};
vm.createContext(ctx); vm.runInContext(script,ctx);
ctx.api=async(path)=>({"/api/conditions":{mv_conc:50,dv_conc:3,dntp_conc:0.8,dna_conc:200},"/api/report":{html:"<html>r</html>",csv:"a,b",n_assays:5},"/api/multiplex":{n_assays:5,n_oligos:14,threshold:-9,channel_conflicts:[{dye:"FAM",assays:["IFNG","IL4"]}],cross_dimers:[],n_flagged:0}}[path]||{});
(async()=>{let ok=0,bad=[];const T=async(l,fn)=>{try{await fn();ok++;console.log("  \u2713 "+l);}catch(e){bad.push(l);console.log("  \u2717 "+l+" -> "+e.message);}};
 await T("initConditions (populate fields)",async()=>{await ctx.initConditions();if(el("cond_dv").value!=3&&el("cond_dv").value!=="3")throw new Error("dv not set");});
 await T("applyConditions (POST + persist + note)",async()=>{await ctx.applyConditions();if(!store["of_cond"])throw new Error("not persisted");});
 await T("seedFSJ (loads 5 locked assays)",async()=>{ctx.seedFSJ();const h=el("wk_out").innerHTML;if(["IFNG","IL4","RPL13","YWHAZ","Plasmodium"].some(g=>h.indexOf(g)<0))throw new Error("missing assay in panel");});
 await T("genReport after seed",async()=>{await ctx.genReport();});
 await T("checkMultiplex after seed (all FAM -> conflict)",async()=>{await ctx.checkMultiplex();if(el("wk_mx").innerHTML.indexOf("conflict")<0)throw new Error("no conflict shown");});
 console.log("\n"+ok+" ok, "+bad.length+" failed"+(bad.length?": "+bad.join(", "):""));
})();
