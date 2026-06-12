import { ParameterPanel } from "./ParameterPanel";
import { MapView } from "./MapView";
import { HazardCurveChart } from "./HazardCurveChart";
import { ResultsTable } from "./ResultsTable";

export function Dashboard() {
  return (
    <div className="h-full flex flex-col xl:flex-row gap-6">
      {/* Left sidebar for parameters */}
      <div className="xl:w-80 shrink-0 flex flex-col gap-4">
        <ParameterPanel />
      </div>

      {/* Main visualization area */}
      <div className="flex-1 flex flex-col gap-6 min-w-0 h-full overflow-hidden">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 h-[400px] shrink-0 min-h-0">
          <MapView />
          <HazardCurveChart />
        </div>
        <div className="flex-1 min-h-[300px] min-w-0">
          <ResultsTable />
        </div>
      </div>
    </div>
  );
}
