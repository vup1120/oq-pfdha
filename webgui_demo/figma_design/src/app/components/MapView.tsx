import { Layers, Maximize2, Plus, Minus } from "lucide-react";
import mapImage from "../../imports/image.png";
import { ImageWithFallback } from "./figma/ImageWithFallback";

export function MapView() {
  return (
    <div className="bg-white rounded-lg border border-slate-200 shadow-sm flex flex-col overflow-hidden relative group">
      <div className="absolute top-4 right-4 z-10 flex flex-col gap-2">
        <div className="bg-white rounded-md shadow-md flex flex-col border border-slate-200">
          <button className="p-2 hover:bg-slate-50 border-b border-slate-100 text-slate-600 transition-colors">
            <Plus className="w-4 h-4" />
          </button>
          <button className="p-2 hover:bg-slate-50 text-slate-600 transition-colors">
            <Minus className="w-4 h-4" />
          </button>
        </div>
        <button className="bg-white rounded-md shadow-md p-2 border border-slate-200 hover:bg-slate-50 text-slate-600 transition-colors">
          <Layers className="w-4 h-4" />
        </button>
      </div>

      <div className="absolute bottom-4 left-4 z-10">
        <div className="bg-white/90 backdrop-blur-sm px-3 py-1.5 rounded-md shadow-md border border-slate-200 text-xs font-medium text-slate-700 flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-red-500 inline-block border border-white shadow-sm"></span>
          Site Location
          <span className="w-4 h-0.5 bg-orange-500 inline-block ml-2"></span>
          Fault Line
        </div>
      </div>

      <div className="absolute top-4 left-4 z-10 opacity-0 group-hover:opacity-100 transition-opacity">
        <button className="bg-white/90 backdrop-blur-sm rounded-md shadow-md p-2 border border-slate-200 hover:bg-slate-50 text-slate-600 transition-colors">
          <Maximize2 className="w-4 h-4" />
        </button>
      </div>

      <div className="flex-1 w-full relative bg-slate-100 overflow-hidden">
        {/* We use ImageWithFallback to render the imported image */}
        <ImageWithFallback 
          src={mapImage} 
          alt="Map View" 
          className="w-full h-full object-cover"
        />
        
        {/* Overlay elements simulating map data */}
        <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2">
          {/* Fault line overlay */}
          <div className="w-48 h-0.5 bg-orange-500 -rotate-45 absolute top-0 left-[-50px] shadow-[0_0_8px_rgba(249,115,22,0.8)]"></div>
          
          {/* Site marker */}
          <div className="absolute top-[20px] left-[30px] flex flex-col items-center pointer-events-none">
            <div className="w-4 h-4 rounded-full bg-red-500 border-2 border-white shadow-md animate-pulse"></div>
            <div className="mt-1 bg-white px-2 py-0.5 rounded text-[10px] font-bold shadow-sm whitespace-nowrap">Site A</div>
          </div>
        </div>
      </div>
    </div>
  );
}
