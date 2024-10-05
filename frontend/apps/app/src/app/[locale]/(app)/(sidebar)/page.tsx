// import { ChartSelectors } from "@/components/charts/chart-selectors";
// import { Charts } from "@/components/charts/charts";
// import { EmptyState } from "@/components/charts/empty-state";
// import { OverviewModal } from "@/components/modals/overview-modal";
// import { Widgets } from "@/components/widgets";
// import { Cookies } from "@/utils/constants";
// import {
//   getBankAccountsCurrencies,
//   getTeamBankAccounts,
// } from "@midday/supabase/cached-queries";
import { cn } from "@v1/ui";
// import { startOfMonth, startOfYear, subMonths } from "date-fns";
import type { Metadata } from "next";
import { cookies } from "next/headers";

export const metadata: Metadata = {
  title: "Overview | V1",
};

// const defaultValue = {
//   from: subMonths(startOfMonth(new Date()), 12).toISOString(),
//   to: new Date().toISOString(),
//   period: "monthly",
// };

export default async function Overview() {
  return (
    <>
      <div>
        <div className="h-[530px] mb-4">
          {/* <ChartSelectors defaultValue={defaultValue} currency={currency} />

          <div className="mt-8 relative">
            {isEmpty && <EmptyState />}

            <div className={cn(isEmpty && "blur-[8px] opacity-20")}>
              <Charts
                value={value}
                defaultValue={defaultValue}
                disabled={isEmpty}
                currency={currency}
                type={chartType}
              />
            </div>
          </div>*/}
        </div>

        {/* <Widgets
          initialPeriod={initialPeriod}
          disabled={isEmpty}
          searchParams={searchParams}
        /> */}
      </div>
      {/* 
      <OverviewModal defaultOpen={isEmpty && !hideConnectFlow} /> */}
    </>
  );
}
