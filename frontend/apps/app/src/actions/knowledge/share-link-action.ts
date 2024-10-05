"use server";

import { authActionClient } from "@/actions/safe-action";
import { dub } from "@/lib/dub";
import { shareLinkSchema } from "./schema";

export const shareLinkAction = authActionClient
  .schema(shareLinkSchema)
  .metadata({
    name: "share-link",
  })
  .action(async ({ parsedInput: { knowledgeId, baseUrl } }) => {
    const link = await dub.links.create({
      url: `${baseUrl}/knowledge/${knowledgeId}`,
    });

    return link?.shortLink;
  });
