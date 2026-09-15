/**
 * Database seeding is disabled per user request to maintain an empty database
 * and maintain only the recent 100-query CSV logbook.
 */
export async function seedDatabase(): Promise<void> {
  // No-op: Seeding disabled to keep database clean.
  return Promise.resolve();
}
