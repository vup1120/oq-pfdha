import { ResponsiveContainer, LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend } from "recharts";
import { Download } from "lucide-react";

const mockData = [
  { displacement: 0.1, petersen: 1e-2, youngs: 1.2e-2, combined: 1.1e-2 },
  { displacement: 0.5, petersen: 5e-3, youngs: 6e-3, combined: 5.5e-3 },
  { displacement: 1.0, petersen: 1e-3, youngs: 1.5e-3, combined: 1.2e-3 },
  { displacement: 2.0, petersen: 5e-4, youngs: 8e-4, combined: 6.5e-4 },
  { displacement: 5.0, petersen: 1e-4, youngs: 2e-4, combined: 1.5e-4 },
  { displacement: 10.0, petersen: 2e-5, youngs: 5e-5, combined: 3.5e-5 },
  { displacement: 20.0, petersen: 5e-6, youngs: 1e-5, combined: 7.5e-6 },
  { displacement: 50.0, petersen: 1e-6, youngs: 2e-6, combined: 1.5e-6 },
];

const formatScientific = (tickItem: number) => {
  if (tickItem === 0) return "0";
  return tickItem.toExponential(1); // Use 1 decimal place to prevent duplicate string keys
};

export function HazardCurveChart() {
  return (
    <div className="bg-white rounded-lg border border-slate-200 shadow-sm flex flex-col overflow-hidden h-full min-w-0">
      <div className="p-4 border-b border-slate-200 flex items-center justify-between shrink-0">
        <h3 className="font-semibold text-slate-800 text-sm">Hazard Curves</h3>
        <button className="text-slate-500 hover:text-slate-800 transition-colors" title="Export Chart Data">
          <Download className="w-4 h-4" />
        </button>
      </div>
      <div className="flex-1 relative min-h-[300px] min-w-0">
        <div className="absolute inset-0 p-4 pb-6">
          <ResponsiveContainer width="100%" height="100%">
            <LineChart
              data={mockData}
              margin={{ top: 10, right: 30, left: 20, bottom: 20 }}
            >
              <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" vertical={false} />
              <XAxis 
                dataKey="displacement" 
                type="number" 
                domain={[0, 'dataMax']}
                tick={{ fontSize: 12, fill: '#64748b' }}
                label={{ value: 'Displacement (cm)', position: 'bottom', offset: 0, fontSize: 12, fill: '#475569' }}
              />
              <YAxis 
                scale="log" 
                domain={[1e-6, 1e-2]} 
                ticks={[1e-6, 1e-5, 1e-4, 1e-3, 1e-2]}
                tickFormatter={formatScientific}
                tick={{ fontSize: 12, fill: '#64748b' }}
                label={{ value: 'Annual Exceedance Freq.', angle: -90, position: 'insideLeft', offset: -10, fontSize: 12, fill: '#475569' }}
              />
              <Tooltip 
                formatter={(value: number) => value.toExponential(2)}
                labelFormatter={(label) => `Displacement: ${label} cm`}
                contentStyle={{ borderRadius: '8px', border: '1px solid #e2e8f0', fontSize: '12px' }}
              />
              <Legend wrapperStyle={{ fontSize: '12px', paddingTop: '20px' }} />
              <Line 
                type="monotone" 
                dataKey="petersen" 
                name="Petersen et al." 
                stroke="#0f2846" 
                strokeWidth={2}
                dot={{ r: 4, fill: '#0f2846' }} 
                activeDot={{ r: 6 }} 
              />
              <Line 
                type="monotone" 
                dataKey="youngs" 
                name="Youngs et al." 
                stroke="#3b82f6" 
                strokeWidth={2}
                dot={{ r: 4, fill: '#3b82f6' }} 
              />
              <Line 
                type="monotone" 
                dataKey="combined" 
                name="Mean Hazard" 
                stroke="#f97316" 
                strokeWidth={2}
                strokeDasharray="5 5"
                dot={false}
              />
            </LineChart>
          </ResponsiveContainer>
        </div>
      </div>
    </div>
  );
}
