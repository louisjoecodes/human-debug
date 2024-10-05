import { getCaseById } from "@v1/supabase/queries";

export default async function Page({ params }: { params: { id: string } }) {
    const { data } = await getCaseById(params.id);

    return (
        <div>
            <h1>{params.id}</h1>
            <pre>{JSON.stringify(data, null, 2)}</pre>
        </div>
    );
}

